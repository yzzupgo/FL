import copy
from fileinput import hook_encoded
import itertools
import json
import math
import os
import random
import subprocess
import time  # 导入 time 模块

def read_json_file(file_path):
    """
    读取 JSON 文件，并返回解析后的 Python 字典。
    
    参数：
    - file_path: JSON 文件的路径
    
    返回：
    解析后的 Python 字典
    """
    try:
        with open(file_path, 'r') as file:
            data = json.load(file)
            return data
    except FileNotFoundError:
        print(f"文件 '{file_path}' 未找到。")
        return None
    except json.JSONDecodeError as e:
        print(f"JSON 解析错误：{e}")
        return None

def generate_combinations(arrays, N):
    # print(f'arrayd {arrays}')
    # print(arrays, faultyLineMutantIndex) 
    # print(faultyLineMutantIndex)
    # sum = []
    # tmp = []
    # for index, i in enumerate(arrays):
    #     for j in faultyLineMutantIndex:
    #         if j in i:
    #             if index not in sum:
    #                 sum.append((index, len(i)))   
    #             break
    # # print(sum)
    # if len(sum) == 0:
    #     return tmp
    # for j in faultyLineMutantIndex:
    #     num = 0
    #     for index, i in enumerate(arrays):
    #         if j in i:
    #             continue
    #         if num >= 66:
    #             break
    #         for k in i:
    #             if num >= 66:
    #                 break
    #             num += 1
    #             tmp.append((j, k))

    # # print(tmp)
    # # print(len(faultyLineMutantIndex))   
    # return tmp
    # 第一步 随机选择两个类簇
    n = len(arrays)
    indices = set(range(n))
    # print(type(indices))
    combinations = []
    result = []
    lastIndex = -1
    lastIndexContent = []

    # 如果数量是奇数，选择一个元素先进行组合
    if n % 2 == 1:
        index1 = random.choice(list(indices))
        indices.remove(index1)
        index2 = random.choice(list(indices))
        lastIndex = index2
        for i in arrays[index2]:
            lastIndexContent.append(i)
        # print(index1, index2)
        combination = (index1, index2)
        combinations.append(combination)
        # print(combination)
    while indices:
        # 随机选择一个索引
        index1 = random.choice(list(indices))
        indices.remove(index1)

        # 随机选择另一个不同的索引
        index2 = random.choice(list(indices))
        indices.remove(index2)

        # 组合对应的元素，并添加到结果中
        combination = (index1, index2)
        combinations.append(combination)
    # print(combinations)
    # 第二步 对于随机选择出的类簇组合再进行组合
    # print(lastIndex, lastIndexContent)
    for i, j in combinations:
        # print(type(i))
        if (i == lastIndex or j == lastIndex) and (len(arrays[i]) == 0 or len(arrays[j]) == 0):
            arrays[lastIndex] = lastIndexContent
        tmp1 = arrays[i] if len(arrays[i]) > len(arrays[j]) else arrays[j]
        tmp2 = arrays[j] if len(arrays[i]) > len(arrays[j]) else arrays[i]
        sub = len(tmp1) - len(tmp2)
        # print(tmp1, tmp2, i, j)

        # 从数组中随机挑选 k 个元素
        selected_elements = random.sample(tmp1, sub)

        # 从原数组中删除这些元素
        tmp1 = [element for element in tmp1 if element not in selected_elements]
        while len(selected_elements) > 0:
            index1 = random.choice(selected_elements)
            selected_elements.remove(index1)
            index2 = random.choice(tmp2)
            result.append((index1, index2))
        while len(tmp1) > 0:
            index1 = random.choice(tmp1)
            tmp1.remove(index1)

            index2 = random.choice(tmp2)
            tmp2.remove(index2)
            result.append((index1, index2))
    # print(result)
    # print(len(result))
    return result
    # Step 1: 计算每个数组的总元素个数
    array_lengths = [len(arr) for arr in arrays]
    total_length = sum(array_lengths)

    # Step 2: 计算每个数组应该抽取的元素个数
    extraction_counts = []
    if N <= len(arrays):
        extraction_counts = [1 for length in array_lengths]
    elif total_length <= N:
        extraction_counts = [length for length in array_lengths]
    else:
        num = N // len(arrays)
        extraction_counts = [min(num, length) for length in array_lengths]
        while sum(extraction_counts) < N:
            sub = N - sum(extraction_counts)
            for i, v in enumerate(extraction_counts):
                if sub < 0:
                    break
                if v < array_lengths[i]:
                    extraction_counts[i] += 1
                    sub -= 1
    # print(extraction_counts, sum(extraction_counts))
    # Step 3: 从每个数组中抽取元素
    selected_elements = []
    for i, arr in enumerate(arrays):
        for item in faultyLineMutantIndex:
            if extraction_counts[i] <= 0:
                break
            if item in arr:
                selected_elements.append(item)
                arr.remove(item)
                extraction_counts[i] -= 1
        selected_elements.extend(random.sample(arr, extraction_counts[i]))

    # Step 4: 使用 itertools.combinations 生成所有可能的两两组合
    combinations = list(itertools.combinations(selected_elements, 2))
    return combinations

def solve_for_n(N):
    discriminant = 1 + 8 * N
    if discriminant < 0:
        return None  # 无实数解

    sqrt_discriminant = math.sqrt(discriminant)
    n1 = (-1 + sqrt_discriminant) / 2
    n2 = (-1 - sqrt_discriminant) / 2

    # 寻找最接近的两个整数解
    closest_n1 = round(n1)
    closest_n2 = round(n2)

    # 选择距离 N 更接近的整数解
    if abs(N - closest_n1) < abs(N - closest_n2):
        return closest_n1
    else:
        return closest_n2

def init(pid, svid, evid):
    # juleiCenterData值为mutantData中的变异体对应的下标
    # juleiCenterPath = f"/home/changzexing/MajorJuleiInfo/{pid}/juleiInfoShanjianban.json"
    # juleiCenterData = read_json_file(juleiCenterPath)
    total_success = 0  # 用于统计所有版本的成功次数
    total_fail = 0     # 用于统计所有版本的失败次数
    total_HOM = 0      # 用于统计生成的高阶变异体总数
    total_time = 0     # 用于统计所有版本生成高阶变异体的总耗时
    for vid in range(svid, evid + 1):
        # 所有的变异体，格式source-org-jfree-chart-renderer-category-AbstractCategoryItemRenderer-235-2
        mutantPath = f"/home/changzexing/MajorJuleiInfo/{pid}/{vid}b/mutantInfo.json"
        if not os.path.exists(mutantPath):
            print(f'{pid}-{vid} 不存在')
            continue
        print(f"{pid}-{vid}")
        success = fail = 0
        HOM_count = 0  # 当前版本尝试生成的高阶变异体数量
        if not os.path.exists(mutantPath):
            continue
        mutantData = read_json_file(mutantPath)
        # if pid != 'Math':
        #     for index, val in enumerate(mutantData):
        #         mutantData[index] = mutantData[index].split('-')
        #         mutantData[index].insert(-1, mutantData[index][-2])
        #         mutantData[index] = '-'.join(mutantData[index])
        # 变异体的详细信息，格式
        #  "source-org-jfree-chart-renderer-category-AbstractCategoryItemRenderer-79-1": {
        #     "index": 1,
        #     "linenum": 79,
        #     "typeOp": "STD",
        #     "mutFilePath": "/home/fanluxi/pmbfl/mutantsFile/Chart/1b/1/org/jfree/chart/renderer/category/AbstractCategoryItemRenderer.java",
        #     "relativePath": "source/org/jfree/chart/renderer/category/AbstractCategoryItemRenderer.java"
        # }
        mutantMapPath = f"/home/changzexing/MaxMajorKillJuzhen/{pid}/{vid}b/mapmajor.json"
        mutantMapData = read_json_file(mutantMapPath)
        if mutantMapData is None or not isinstance(mutantMapData, dict):
            print(f"mutantMapData 为空或不是字典，跳过版本 {pid}-{vid}")
            continue
        if mutantMapData is None:
            mutantMapData = []
        num = len(mutantData)
        # faultyLinePath = f"/home/changzexing/faultyline/faultyLine/{pid}/{pid}_{vid}_FalutLine.json"
        # if not os.path.exists(faultyLinePath):
        #     continue
        # faultyLineData = read_json_file(faultyLinePath)
        # faultyLineMutant = []
        # faultyLineMutantIndex = []
        # for item in faultyLineData:
        #     if len(faultyLineData[item]) == 0:
        #         continue
        #     for _item in faultyLineData[item]:
        #         faultyLineMutant.append(f"{item[1:-5].replace('/', '-')}-{_item}")
        # # print(faultyLineMutant)
        # for index, value in enumerate(mutantData):
        #     tmp = '-'.join(value.split('-')[:-1])
        #     if tmp in faultyLineMutant:
        #         faultyLineMutantIndex.append(index)
                # print(index, mutantData[index], tmp)
        # print(faultyLineMutantIndex, mutantData[faultyLineMutantIndex[0]])
        # break
        juleiCenterPath = f"/home/changzexing/MajorJuleiInfo/{pid}/{vid}b/juleiInfo.json"
        juleiCenterData = read_json_file(juleiCenterPath)
        vidJuleiCenter = juleiCenterData
        # print(len(mutantData), len(vidJuleiCenter))
        data = []
        for item in vidJuleiCenter:
                data.append(list(vidJuleiCenter[item]))
        # print(data)
        # print(vidJuleiCenter)
        # selectNum = solve_for_n(num)
        if len(data) < 2:
            continue
        HOMInfo = []
        HOMInfo = generate_combinations(data, num)
        HOM_count = len(HOMInfo)  # 记录当前版本的高阶变异体组合数量
        total_HOM += HOM_count    # 累加到总的高阶变异体数量
        print(f"{pid}-{vid} 生成的高阶变异体组合数量：{HOM_count}")
        # print(num, len(HOMInfo))
        # print(HOMInfo) 
        # continue
        print(len(HOMInfo))
        # continue
        flagHOM = {}
        start_time = time.time()  # 开始计时
        for m1, m2 in HOMInfo:
            # print(mutantData[m1], mutantData[m2])
            # print(mutantMapData[mutantData[m1]], mutantMapData[mutantData[m1]])
            tmp = mutantData[m1].split('-')
            # if pid != 'Math' and tmp[-2] == tmp[-3]:
            #     filem1 = '/'.join(mutantData[m1].split('-')[:-3])
            #     filem2 = '/'.join(mutantData[m1].split('-')[:-3])
            #     linem1 = mutantData[m1].split('-')[-3]
            #     linem2 = mutantData[m2].split('-')[-3]
            # else:
            filem1 = '/'.join(mutantData[m1].split('-')[:-2])
            filem2 = '/'.join(mutantData[m1].split('-')[:-2])
            linem1 = mutantData[m1].split('-')[-2]
            linem2 = mutantData[m2].split('-')[-2]
            # print(mutantData[m1], mutantData[m2])
            if filem1 != filem2:
                continue
            if filem1 not in flagHOM:
                flagHOM[filem1] = {}
            if f"{linem1}+{linem2}" not in flagHOM[filem1]:
                flagHOM[filem1][f"{linem1}+{linem2}"] = 1
            else:
                flagHOM[filem1][f"{linem1}+{linem2}"] += 1
            # if pid != 'Math':
            #     filename = mutantData[m1].split('-')[-4]
            # else:
            filename = mutantData[m1].split('-')[-3]
            # print(filename)
            # print(filem1, filem2)
            # print(flagHOM)
            zuhepath = f"/home/changzexing/zuhemutantMajorSFClu/{pid}/{pid.lower()}_{vid}_buggy/{filem1}/{linem1}+{linem2}/{flagHOM[filem1][f'{linem1}+{linem2}']}"
            # print(zuhepath)
            if not os.path.exists(zuhepath):
                os.makedirs(zuhepath)
            zuhepath += f"/{filename}.java"
            if pid == 'Math':
                path1 = mutantMapData[mutantData[m1]]["mutFilePath"].replace("/home/fanluxi/pmbfl/mutantsFile", "/home/changzexing/mutantFaultyFileMajor")
            else:
                path1 = mutantMapData[mutantData[m1]]["mutFilePath"]
                # path1 = mutantMapData[mutantData[m1]]["mutFilePath"].replace('mutantFautyfileMajor', 'mutantFaultyFileMajor')
            if pid == 'Math':
                path2 = mutantMapData[mutantData[m2]]["mutFilePath"].replace("/home/fanluxi/pmbfl/mutantsFile", "/home/changzexing/mutantFaultyFileMajor")
            else:
                path2 = mutantMapData[mutantData[m2]]["mutFilePath"]
                # path2 = mutantMapData[mutantData[m2]]["mutFilePath"].replace('mutantFautyfileMajor', 'mutantFaultyFileMajor')
            originpath = f"/home/changzexing/d4jclean/{pid}/{vid}b/{filem1}.java"
            # print(vj[0], vj[1], v, vl)
            # print(path1, path2)
            # print(zuhepath)
            # print(originpath, os.path.exists(originpath))
            # print(path1, os.path.exists(path1))
            # print(path2, os.path.exists(path2))
            # break
            command = f"""
java -cp /home/changzexing/d4jscript/codeMerge-1.0-SNAPSHOT.jar:/home/changzexing/d4jscript/javaparser-core-3.24.0.jar:/home/changzexing/d4jscript/java-diff-utils-4.12.jar:/home/changzexing/d4jscript/commons-cli-1.4.jar org.example.Main -pathOri {originpath} -pathM1 {path1} -pathM2 {path2} -pathOutput {zuhepath}
"""
            copyCodeFlag = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if copyCodeFlag.returncode == 0:
                print(f"{zuhepath}成功")
                # print(copyCodeFlag.stdout)
                success += 1
            else:
                print(f"{zuhepath}失败")
                print(copyCodeFlag.stderr)
                fail += 1
            # break
        # break
        end_time = time.time()    # 结束计时
        time_spent = end_time - start_time
        total_time += time_spent  # 累加总耗时
        print(f"{pid}-{vid} 完成，尝试生成高阶变异体数量：{HOM_count}，成功：{success}，失败：{fail}")
        total_success += success
        total_fail += fail
        print(f"{pid}-{vid} 完成 成功{success} 失败{fail}")
    print(f"项目 {pid} 总共尝试生成高阶变异体数量：{total_HOM}，成功：{total_success}，失败：{total_fail}，总耗时：{total_time:.2f}秒")

if __name__ == '__main__':
    # init("Chart", 1, 26)
    # init('Time', 1, 26)
    # init('Math', 1, 1)
    # init('Lang', 1, 8)
    # init('Math', 1, 106)
    init('Closure', 1, 5)
    # init('Cli', 4, 4)
    # init('Codec', 1, 18)
    # init('Csv', 1, 16)
    # init('Mockito', 1, 38)