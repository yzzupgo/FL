from concurrent.futures import ProcessPoolExecutor
import Last2First
import RandomMix
import RandomMixMajorMbert
import Last2FirstMajorMbert

def init(pid, len, index, processName):
    mp = {
        'Chart': 26,
        'Time': 27,
        'Lang': 65,
        'Math': 106,
        'Mockito': 38
    }
    if index >= len:
        # RandomMixMajorMbert.init(pid, 1, mp[pid], 'RandomMixMajorMbert', 'tmpd4jclean')
        RandomMix.init(pid, 1, mp[pid])
        # print(pid, 1, mp[pid], 'RandomMixMajorMbert', 'tmpd4jclean', processName)
    else:
        # Last2FirstMajorMbert.init(pid, 1, mp[pid], 'Last2FirstMajorMbert', 'd4jclean')
        Last2First.init(pid, 1, mp[pid])
        # print(pid, 1, mp[pid], 'Last2FirstMajorMbert', 'd4jclean', processName)

if __name__ == '__main__':
    # Last2First.init('Math', 1, 106)
    # RandomMix.init('Math', 1, 106)
    # versionList = ['Chart', 'Time', 'Lang']
    versionList = ['Mockito']
    len1 = len(versionList)
    num = 2 * len1
    with ProcessPoolExecutor(max_workers=num) as pool:
        # Submit the jobs to the pool
        futures = [pool.submit(init, versionList[i % len1], len1, i, 'process%d'%(i,)) for i in range(num)]
        
        # Wait for all the jobs to complete
        for future in futures:
            future.result()