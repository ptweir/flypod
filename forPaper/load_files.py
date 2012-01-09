import os, pickle

def compare_dir_names(dn1, dn2):
    n1, n2 = int(dn1[3:]),int(dn2[3:])
    return cmp(n1,n2)
def compare_file_times(fn1, fn2):
    t1, t2 = int(fn1[4:12]+fn1[13:19]),int(fn2[4:12]+fn2[13:19])
    return cmp(t1,t2)
def load_flies_skies(rootDir, baseDirs):
    skies, flies = {}, {}
    for b, baseDir in enumerate(baseDirs):
        sks, fls = {}, {}
        dNames = os.listdir(os.path.join(rootDir,baseDir))
        flyDirNames = [D for D in dNames if D[:3] == 'fly']
        flyDirNames.sort(compare_dir_names)

        for dNum, dName in enumerate(flyDirNames):
            dirName=os.path.join(rootDir,baseDir,dName)
            print dirName
            
            filenames = os.listdir(dirName)
            pklFlyFilenames = [f for f in filenames if f[:4] == 'aFly' and f[-3:] == 'pkl']
            pklFlyFilenames.sort(compare_file_times)
            pklFlyFilename = pklFlyFilenames[-1]
            inPklFlyFile = open(os.path.join(dirName,pklFlyFilename), 'rb')
            fly = pickle.load(inPklFlyFile)
            inPklFlyFile.close()
            
            pklSkyFilenames = [f for f in filenames if f[:4] == 'aSky' and f[-3:] == 'pkl']
            pklSkyFilenames.sort(compare_file_times)
            pklSkyFilename = pklSkyFilenames[-1]
            inPklSkyFile = open(os.path.join(dirName,pklSkyFilename), 'rb')
            sky = pickle.load(inPklSkyFile)
            inPklSkyFile.close()
            
            fls[dNum] = fly.copy()
            sks[dNum] = sky.copy()
            
        flies[baseDir] = fls
        skies[baseDir] = sks
    return flies, skies

rootDir = '/home/peter/data/'

baseDirs = ['indoor/dark','indoor/halogen','blueFilter','circularPolarizer', 'grayFilter', 'noFilter']

flies, skies = load_flies_skies(rootDir, baseDirs)
