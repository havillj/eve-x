from util import *

def main():
    start = 1
    end = 100000
    if len(sys.argv) >= 4:
        start = int(sys.argv[3])
    if len(sys.argv) >= 5:
        end = int(sys.argv[4])
        
    print(getContig(SPECIMENS_DIR + '/' + sys.argv[1], sys.argv[2], start, end))

main()
