from pathlib import Path
import os

# def fuck():
#     dir = Path('/Volumes/Data/aegypti/all/').resolve()
#     for subdir in dir.iterdir():
#         if not subdir.is_dir():
#             print(subdir)
#             os.system('mv ' + str(subdir) + ' ' + str(subdir) + '.tmp')
#             os.mkdir(str(subdir))
#             os.system('mv ' + str(subdir) + '.tmp ' + str(subdir) + '/' + subdir.name + '.sorted.deduped.merged.bam')
#
#             print('mv ' + str(subdir) + ' ' + str(subdir) + '.tmp')
#             print('mkdir ' + str(subdir))
#             print('mv ' + str(subdir) + '.tmp ' + str(subdir) + '/' + subdir.name + '.sorted.deduped.merged.bam')
#             print()

def main():
    dir = Path('/Volumes/Data/aegypti/all/').resolve()
    for subdir in dir.iterdir():
        if subdir.is_dir():
            try:
                subdir.name.index('.')  # if subdir name contains '.' then skip
                continue
            except:
                for file in subdir.iterdir():
                    try:
                        index = file.name.index('.sorted.deduped.merged.bam')
                    except:
                        print('file name: ' + file.name)
                        continue
                    try:
                        index2 = file.name.index('-aegy-')
                        index2 += 6
                    except:
                        try:
                            index2 = file.name.index('_aegypti_')
                            index2 += 9
                        except:
                            try:
                                index2 = file.name.index('-masc-')
                                index2 += 6
                            except:
                                index2 = 0
                    newDir = file.name[index2:index]
                    if not Path('/Volumes/Data/aegypti/analyzed/' + newDir).exists():
                        print(newDir)
                        try:
                            os.mkdir('/Volumes/Data/aegypti/all/' + newDir)
                        except:
                            pass
                        os.system('mv '+ str(file.resolve()) + ' /Volumes/Data/aegypti/all/' + newDir)

main()
