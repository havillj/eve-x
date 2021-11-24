# kmeans.py

from util import *
    
def distance(p, q):  # identities
    return np.count_nonzero(p - q)  # same as len(p) - np.count_nonzero((p - q) == 0)
    
# def getFurthestPoint(centroids, data):
#     """ Return random point with probabilities proportional to
#         minimum distance from all current centroids."""
#     
#     if len(centroids) == 0:
#         return data[random.randrange(len(data))]
#         
#     prob = np.zeros(len(data))
#     sumDist = 0
#     for index in range(len(data)):
#         minDist = distance(data[index], centroids[0])
#         for centroid in centroids[1:]:
#             minDist = min(minDist, distance(data[index], centroid))
#         sumDist += minDist
#         prob[index] = sumDist
#             
#     prob /= sumDist
#         
#     r = random.random()
#     for index in range(len(data)):
#         if r < prob[index]:
#             return data[index]
    
def centroid(cluster, data, centroids):
    n = len(cluster)
    if n == 0:
        theCentroid = data[random.randrange(len(data))]
        while len(centroids) < len(np.unique(data, axis=0)) and any((theCentroid == c).all() for c in centroids):
            theCentroid = data[random.randrange(len(data))]
        return theCentroid
#      return getFurthestPoint(centroids, data)
    if n == 1:
        return data[cluster[0]]
    m = len(data[0])   # dimension of points
    clusterData = [data[index] for index in range(len(data)) if index in cluster]
    theCentroid = np.round(np.sum(clusterData, axis=0) / n)
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
            
        if sumDist > 0:
            prob /= sumDist
            
            r = random.random()
            for index in range(len(data)):
                if r < prob[index]:
                    centroidIndices.append(index)
                    centroids.append(data[index])
        else:  # remaining points are identical to an existing centroid, so just pick one
            for index in range(len(data)):
                if index not in centroidIndices:
                    centroidIndices.append(index)
                    centroids.append(data[index])
                    break
    
    return centroids

def kmeans(data, k):
    n = len(data)
    centroids = initialCentroids(data, k)  # kmeans++
    
    prevClusters = None
    clusters = []
    iteration = 0
    while (iteration < config['MAX_KMEANS_ITERATIONS']) and (clusters != prevClusters):
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
        
    return clusters
    
def similarity(pIndex, cluster, data):
    return sum([distance(data[pIndex], data[index]) for index in cluster if index != pIndex]) / (len(cluster) - 1)
    
def dissimilarity(pIndex, clusters, data):
    minDissimilarity = len(data[0])
    for cluster in clusters:
        if len(cluster) > 0 and pIndex not in cluster:
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

def findClusters(fileName, omitFirst = True, minK = config['MIN_K'], maxK = config['MAX_K']):
    print('Clustering aligned sequences in ' + fileName + '...')
    
    aln = AlignIO.read(fileName, 'fasta') # Bio.Align.MultipleSeqAlignment object
    
    if len(aln) <= 2:
        print('  Skipping - too few sequences to cluster.')
        return None

    d = {'-': 0, 'a': 1, 'c': 2, 'g': 3, 't': 4, 'A': 1, 'C': 2, 'G': 3, 'T': 4}
    
    refSeq = aln[0].seq
    data = []
    ids = []
    first = len(refSeq) - 1
    last = 0
    for i in range(int(omitFirst), len(aln)):
        seq = aln[i].seq
        vector = np.zeros(len(refSeq))
        for j in range(len(seq)):
            if seq[j] != '-':
                first = min(first, j)
                last = max(last, j)
            vector[j] = d.get(seq[j], 5)
        data.append(vector)
        ids.append(aln[i].id)
    
    for i in range(len(data)):
        data[i] = data[i][first:last+1]

    silhouettes = {}
    clusters = [list(range(len(data)))]
#    silhouettes[1] = [(meanSilhouetteValue(clusters, data), clusters)]
    for k in range(minK, min(maxK, len(data)) + 1):
        silhouettes[k] = []
        print('  k = ' + str(k))
        print('  Trial ', end='')
        for trial in range(1, config['KMEANS_TRIALS'] + 1):
            print(str(trial), end=' ', flush=True)
            clusters = kmeans(data, k)
            silhouettes[k].append((meanSilhouetteValue(clusters, data), clusters))
        print()

    print('\n  Silhouette values:')
    maxSV = -1
#    maxSV = silhouettes[1][0][0]
#    bestK = 1
#    clusters = silhouettes[1][0][1]
    for k in range(minK, min(maxK, len(data)) + 1):
        silhouettes[k].sort(reverse = True)
        print('  {0:>2}: '.format(k), end='')
        sumSV = 0
        for i in range(config['KMEANS_TRIALS'] - 1):
            print('{0:6.4f}, '.format(silhouettes[k][i][0]), end='')
            sumSV += silhouettes[k][i][0]
        sumSV += silhouettes[k][-1][0]
        print('{0:6.4f} (average = {1:6.4f})'.format(silhouettes[k][-1][0], sumSV / config['KMEANS_TRIALS']))
        if silhouettes[k][0][0] > maxSV:
            maxSV = silhouettes[k][0][0]
            bestK = k
            clusters = silhouettes[k][0][1]
            
    realK = sum(len(cluster) > 0 for cluster in clusters)
            
    if realK == bestK:
        print('\n  Best k = ' + str(bestK) + ' clusters have silhouette value = ' + '{0:6.4f}'.format(maxSV) + '.\n')
    else:
        print('\n  Best k = ' + str(realK) + ' (' + str(bestK) + ' minus ' + str(bestK - realK) + ' empty) clusters have silhouette value = ' + '{0:6.4f}'.format(maxSV) + '.\n')

    records = [aln[0]]
    clusters.sort(key = lambda c: len(c), reverse = True)  # start with largest cluster
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

#        print('\n  Best k = ' + str(k) + ' clusters:\n')

    clusteredFileName = fileName.split('.fasta')[0] + '_clustered_k' + str(realK)
    clusterNum = 1
    regionCounts = {}   
    outputText = open(clusteredFileName + '.txt', 'w')
    for cluster in clusters:
        if len(cluster) > 0:
            regionCounts[clusterNum] = {}
            cluster.sort(key = lambda i: ids[i])
            print('  Cluster ' + str(clusterNum))
            outputText.write('  Cluster ' + str(clusterNum) + '\n')
            for index in cluster:
                print('    ' + ids[index])
                outputText.write('    ' + ids[index] + '\n')
                specimenName = ids[index].split('_|_')[0]
                region, num = getSpecimenLabel(specimenName)
                if region not in regionCounts[clusterNum]:
                    regionCounts[clusterNum][region] = [specimenName]
                else:
                    regionCounts[clusterNum][region].append(specimenName)
                if omitFirst:
                    records.append(copy.deepcopy(aln[index + 1]))
                else:
                    records.append(copy.deepcopy(aln[index]))
                d = records[-1].description
                records[-1].description = 'C' + str(clusterNum) + '_|_' + d
                records[-1].id = records[-1].description
                records[-1].name = records[-1].description
            print('\n    Region counts: ')
            outputText.write('\n    Region counts (total / unique specimens): \n')
            sortedRegions = list(regionCounts[clusterNum].keys())
            sortedRegions.sort()
            for region in sortedRegions:
                print('      ' + region + ' - ' + str(len(regionCounts[clusterNum][region])) + ' / ' + str(len(set(regionCounts[clusterNum][region]))))
                outputText.write('      ' + region + ' - ' + str(len(regionCounts[clusterNum][region])) + ' / ' + str(len(set(regionCounts[clusterNum][region]))) + '\n')
            clusterNum += 1
            print()
            outputText.write('\n')
    outputText.close()
    
    newAlignment = MultipleSeqAlignment(records)
    AlignIO.write(newAlignment, clusteredFileName + '.fasta', 'fasta')
    print('  Clustered alignment written to ' + clusteredFileName + '.fasta')
    
    return clusteredFileName + '.fasta'

def usage():
    print('Usage: ' + sys.argv[0] + ' [--min_k=k] [--max_k=k] [--omit_first] filename.fasta')
    
def main():
    if len(sys.argv) < 2:
        usage()
        return
    
    omitFirst = False
    minK = config['MIN_K']
    maxK = config['MAX_K']
    for arg in sys.argv[1:-1]:
        if arg.startswith('--min_k'):
            try:
                minK = int(arg.split('=')[1])
            except:
                usage()
                return
        elif arg.startswith('--max_k'):
            try:
                maxK = int(arg.split('=')[1])
            except:
                usage()
                return
        elif arg == '--omit_first':
            omitFirst = True
        else:
            usage()
            return
    
    if sys.argv[-1][0] == '-':
        usage()
        return
        
    findClusters(sys.argv[-1], omitFirst, minK, maxK)
    
if __name__ == '__main__':
    main()