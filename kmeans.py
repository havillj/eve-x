import math
import random
import copy
import numpy as np
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

MIN_K = 2
MAX_K = 15
KMEANS_TRIALS = 10
MAX_KMEANS_ITERATIONS = 100

# def distance(p, q):
#     sumSquares = 0
#     for i in range(len(p)):
#         sumSquares += (p[i] - q[i]) ** 2
# #    return math.sqrt(sumSquares)
#     return sumSquares

# def distance(p, q):
#     return sum(abs(p - q))
    
# def identities(p, q):
#     count = 0
#     for index in range(len(p)):
#         if p[index] == q[index]:
#             count += 1
#     return -count
    
def distance(p, q):  # identities
#    return len(p) - np.count_nonzero((p - q) == 0)
    return np.count_nonzero(p - q)  # same as above!

# fast way to make this a low distance?
# ------------AAAA-----
# ----XXXXXXXXAAAAXX---

#def distance(p, q):  # identities
#    matchScore = 2
    
#    score = 3 * base matches + 1 * gap matches = 1 * all matches + 2 * base matches
#    3 * len(p) - score = 3 * len(p) - (1 * all matches + 2 * base matches) 
#                       = 3 * (len(p) - all matches) + 2 * all matches - 2 * base matches
#                       = 3 * (len(p) - all matches) + 2 * (all matches - base matches) 
#                       = 3 * (len(p) - all matches) + 2 * (gap matches)
#    len(p) - all matches = np.count_nonzero(p - q)
#    gap matches = np.count_nonzero((p + q) == 0)

#    return matchScore * np.count_nonzero(p - q) + (matchScore - 1) * np.count_nonzero((p + q) == 0)
#    return matchScore * np.count_nonzero(p - q) + (matchScore - 1) * np.count_nonzero(((p + q) > 0) and ((p + q) == p or (p + q) == q)) # syntax
    
def getFurthestPoint(centroids, data):
    if len(centroids) == 0:
        return data[random.randrange(len(data))]
        
    prob = np.zeros(len(data))
    sumDist = 0
    for index in range(len(data)):
        minDist = distance(data[index], centroids[0])
        for centroid in centroids[1:]:
            minDist = min(minDist, distance(data[index], centroid))
        sumDist += minDist
        prob[index] = sumDist
            
    prob /= sumDist
        
    r = random.random()
    for index in range(len(data)):
        if r < prob[index]:
            return data[index]
    
def centroid(cluster, data, centroids):
    n = len(cluster)
    if n == 0:
        return data[random.randrange(len(data))]
#      return getFurthestPoint(centroids, data)
    if n == 1:
        return data[cluster[0]]
    m = len(data[0])   # dimension of points
    clusterData = [data[index] for index in range(len(data)) if index in cluster]
    theCentroid = np.sum(clusterData, axis=0) / n
    return theCentroid
    
def initialCentroids(data, k):
    centroidIndices = [random.randrange(len(data))]
    centroids = [data[centroidIndices[0]]]
    
    while len(centroidIndices) < k:
        prob = np.zeros(len(data))
        sumDist = 0.0
        for index in range(len(data)):
            if index in centroidIndices:
                prob[index] = sumDist
            else:
                minDist = distance(data[index], centroids[0])
                for centroid in centroids[1:]:
                    minDist = min(minDist, distance(data[index], centroid))
                sumDist += minDist
                prob[index] = sumDist
            
        prob /= sumDist
            
        r = random.random()
        for index in range(len(data)):
            if r < prob[index]:
                centroidIndices.append(index)
                centroids.append(data[index])
    
    return centroids

def kmeans(data, k):
    n = len(data)
    centroids = initialCentroids(data, k)  # kmeans++
    
    prevClusters = None
    clusters = []
    iteration = 0
    while (iteration < MAX_KMEANS_ITERATIONS) and (clusters != prevClusters):
        prevClusters = clusters[:]
        clusters = []
        for clustIndex in range(k):
            clusters.append([])
        
        for dataIndex in range(n):
            minIndex = 0
            minDist = distance(centroids[0], data[dataIndex])
            for clustIndex in range(1, k):
                dist = distance(centroids[clustIndex], data[dataIndex])
                if dist < minDist:
                    minIndex = clustIndex
                    minDist = dist
            clusters[minIndex].append(dataIndex)
        
        for clustIndex in range(k):
            centroids[clustIndex] = centroid(clusters[clustIndex], data, centroids[:clustIndex])
        iteration += 1
#         print([len(c) for c in clusters])
#         print([len(c) for c in prevClusters])
        
    return clusters
    
def similarity(pIndex, cluster, data):
#    dataN = np.array([data[index] for index in cluster])
    
#     print('sim')
#     print(sum([distance(data[pIndex], data[index]) for index in cluster if index != pIndex]) / (len(cluster) - 1))
#     print(np.sum(abs(dataN - data[pIndex])) / (len(cluster) - 1))
    
#    return np.sum(abs(dataN - data[pIndex])) / (len(cluster) - 1)

    return sum([distance(data[pIndex], data[index]) for index in cluster if index != pIndex]) / (len(cluster) - 1)
    
def dissimilarity(pIndex, clusters, data):
    minDissimilarity = len(data[0])
    for cluster in clusters:
        if len(cluster) > 0 and pIndex not in cluster:
#            dataN = np.array([data[index] for index in cluster])
#            minDissimilarity = min(minDissimilarity, np.sum(abs(dataN - data[pIndex])) / len(cluster))
            minDissimilarity = min(minDissimilarity, sum([distance(data[pIndex], data[index]) for index in cluster]) / len(cluster))
        
    return minDissimilarity
    
def meanSilhouetteValue(clusters, data):
    total = 0
    for cluster in clusters:
        if len(cluster) > 1:
            for index in cluster:    
                a = similarity(index, cluster, data)
                b = dissimilarity(index, clusters, data)
                total += (b - a) / max(a, b)
    return total / len(data)

def findClusters(fileName):
    aln = AlignIO.read(fileName, 'fasta') # Bio.Align.MultipleSeqAlignment object

    d = {'-': 0, 'a': 1, 'c': 2, 'g': 3, 't': 4, 'A': 1, 'C': 2, 'G': 3, 'T': 4}
    
    refSeq = aln[0].seq
    data = []
    ids = []
    first = len(refSeq) - 1
    last = 0
    for i in range(1, len(aln)):
        seq = aln[i].seq
#         vector = np.zeros(len(refSeq))
#         for j in range(len(seq)):
#             if seq[j] != '-':
#                 first = min(first, j)
#                 last = max(last, j)
#             if seq[j] == refSeq[j]:
#                 vector[j] = 1
# #             elif seq[j] == '-':
# #                 vector[j] = -1
        vector = np.zeros(len(refSeq))
        for j in range(len(seq)):
            if seq[j] != '-':
                first = min(first, j)
                last = max(last, j)
            vector[j] = d.get(seq[j], 6)

        data.append(vector)
        ids.append(aln[i].id)
    
    for i in range(len(data)):
        data[i] = data[i][first:last+1]

    silhouettes = []
    for trial in range(KMEANS_TRIALS):
        print('Trial ' + str(trial))
#        clustersDict = {}
        print('k = ', end='')
        for k in range(MIN_K, MAX_K + 1):
            print(str(k), end=' ', flush=True)
            clusters = kmeans(data, k)
            silhouettes.append((meanSilhouetteValue(clusters, data), k, clusters))
        print()

#     kFreq = {}
#     for sv, k, clusters in bestKs:
#         if k not in kFreq:
#             kFreq[k] = []
#         kFreq[k].append(sv)
#         
#     items = [(len(kFreq[k]), k) for k in kFreq]
#     items.sort(reverse = True)
#     for freq, k in items:
#         avg = sum(kFreq[k]) / len(kFreq[k])
#         print(str(k) + ': ' + str(kFreq[k])[1:-1] + ' (' + str(avg) + ')')

    print(silhouettes)
        
    sv, bestK, clusters = max(silhouettes)

    print('Best k = ' + str(bestK) + ' clusters have silhouette value = ' + str(sv))
    bestClusters = {bestK: clusters}
    triples = [(sv, k, clusters) for (sv, k, clusters) in silhouettes if k == bestK - 1]
    if len(triples) >= 1:
        bestClusters[bestK - 1] = max(triples)[2]
        print(' Best k-1 = ' + str(bestK - 1) + ' clusters have silouette value = ' + str(max(triples)[0]))
    triples = [(sv, k, clusters) for (sv, k, clusters) in silhouettes if k == bestK + 1]
    if len(triples) >= 1:
        bestClusters[bestK + 1] = max(triples)[2]
        print(' Best k+1 = ' + str(bestK + 1) + ' clusters have silouette value = ' + str(max(triples)[0]))
    
    for k in [bestK, bestK - 1, bestK + 1]:
        if k not in bestClusters:
            continue
        clusters = bestClusters[k]
        records = [aln[0]]
        clusters.sort(key = lambda c: len(c), reverse = True)  # start with largest cluster, empty clusters now at end
        centroids = [centroid(cluster, data, []) for cluster in clusters if len(cluster) > 0]
        for index1 in range(1, len(centroids) - 1):
            minDist = distance(centroids[index1 - 1], centroids[index1])
            minIndex = index1
            for index2 in range(index1 + 1, len(centroids)):
                d = distance(centroids[index1 - 1], centroids[index2])
                if d < minDist:
                    minDist = d
                    minIndex = index2
            clusters[index1], clusters[minIndex] = clusters[minIndex], clusters[index1]
            centroids[index1], centroids[minIndex] = centroids[minIndex], centroids[index1]

        print('\nBest k = ' + str(k) + ' clusters:\n')

        clusterNum = 1
        for cluster in clusters:
            cluster.sort(key = lambda i: ids[i])
            print('Cluster ' + str(clusterNum))
            for index in cluster:
                print('  ' + ids[index])
                records.append(copy.deepcopy(aln[index + 1]))
                d = records[-1].description
                records[-1].description = 'C' + str(clusterNum) + '_' + d
                records[-1].id = records[-1].description
                records[-1].name = records[-1].description
            clusterNum += 1
            print()
        newAlignment = MultipleSeqAlignment(records)
        AlignIO.write(newAlignment, fileName.split('.fasta')[0] + '_clustered_k' + str(k) + '.fasta', 'fasta')

#fileName = '/Volumes/Data2/results4/viruses/sequences/Phenuiviridae/sequences_NC_038263.1_per_contig.fasta'
#fileName = '/Volumes/Data2/results4/viruses/sequences/Flaviviridae/sequences_NC_001564.2_per_contig_aligned.fasta'
findClusters('/Volumes/Data2/results4/viruses/sequences/Phenuiviridae/sequences_NC_038263.1_per_contig_aligned.fasta')
#findClusters('/Volumes/Data2/results4/viruses/sequences/Flaviviridae/sequences_NC_001564.2_per_contig_aligned.fasta')
#from msa import drawMSA
#drawMSA(fileName.split('.fasta')[0] + '_clustered.fasta', 'NC_038263.1', '')
