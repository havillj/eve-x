class Hit:
    def __init__(self, pos = [], overlap = False, primary = True, data = None):
        self._pos = pos
        self._overlap = overlap
        self._primary = primary
        self._data = data
        
    def doesOverlap(self):
        return self._overlap
        
    def isPrimary(self):
        return self._primary
        
    def getData(self):
        return self._data
        
    def getPos(self):
        return self._pos
        
    def __len__(self):
        return len(self._pos)
        
    def __getitem__(self, index):
        return self._pos[index]
        
    def __contains__(self, interval):
        try:
            iterator = iter(interval)
        except TypeError:
            start = end = interval
        else:
            start, end = interval
            if start > end:
                start, end = end, start
                
        for index in range(0, len(self._pos), 2):
            if (self._pos[index] <= self._pos[index + 1]) and ((start >= self._pos[index]) and (end <= self._pos[index + 1])):
                return True
            elif (self._pos[index] > self._pos[index + 1]) and ((start <= self._pos[index]) and (end >= self._pos[index + 1])):
                return True
        return False
        
class Hits:
    def __init__(self):
        self._hits = []   # list of hits
        
    def addHit(self, pos, overlap = False, primary = True, data = None):
        self._hits.append(Hit(pos, overlap, primary, data))
        
    def addHits(self, otherHitsList):
        self._hits.extend(otherHitsList)
        
    def getHits(self):
        return self._hits
        
    def getPrimaryHitsOnly(self):
        primaryHits = []
        for hit in self._hits:
            if hit.isPrimary():
                primaryHits.append(hit)
        return primaryHits
        
    def adjust(self, adjustment):
        for index in range(len(self._hits)):
            for index2 in range(len(self._hits[index]._pos)):
                self._hits[index]._pos[index2] += adjustment
                
    def __contains__(self, interval):
        for hit in self._hits:
            if interval in hit:
                return True
        return False
                
    def __len__(self):
        return len(self._hits)
