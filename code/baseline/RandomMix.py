import json
import os
import random
import subprocess

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

def init(pid, svid, evid):
    No = []
    for vid in range(svid, evid + 1):
        if vid in No:
            print(f"{pid}-{vid} 跳过")
            continue
        print(f"{pid}-{vid}")
        success = fail = 0
        mutantInfoPath = f"/home/changzexing/mutant_result_faulty_file_major/{pid}/{vid}b/mutantInfo.json"
        # mutantInfoPath = f"/home/changzexing/mutant_result_faulty_file_major/{pid}/{vid}b/muInfo.json"
        mutantInfoData = read_json_file(mutantInfoPath)
        if mutantInfoData is None:
            print(f"{pid} - {vid} 没有变异体信息")
            continue
        fileMutantInfo = {}
        for index, val in enumerate(mutantInfoData):
            if val['relativePath'] in fileMutantInfo:
                fileMutantInfo[val['relativePath']].append(index)
            else:
                fileMutantInfo[val['relativePath']] = [index]
        num = len(mutantInfoData)
        # print(mutantInfoData[0])
        HOMInfo = []
        for item in fileMutantInfo:
          tmpdata = fileMutantInfo[item]
          while len(tmpdata) >= 2:
            index1 = random.choice(tmpdata)
            tmpdata.remove(index1)
            index2 = random.choice(tmpdata)
            tmpdata.remove(index2)
            HOMInfo.append((index1, index2))
          # print(HOMInfo)  
          print(len(HOMInfo))
        # return
        flagHOM = {}

        for indexm1, indexm2 in HOMInfo:
            m1 = mutantInfoData[indexm1]
            m2 = mutantInfoData[indexm2]
            m1['relativePath'] = '/'.join(m1['relativePath'].split('/')[6:]) 
            filename = m1['relativePath'].split('/')[-1]
            linem1 = m1['linenum']
            linem2 = m2['linenum']
            if m1['relativePath'] not in flagHOM:
                flagHOM[m1['relativePath']] = {}
            if f"{linem1}+{linem2}" not in flagHOM[m1['relativePath']]:
                flagHOM[m1['relativePath']][f"{linem1}+{linem2}"] = 1
            else:
                flagHOM[m1['relativePath']][f"{linem1}+{linem2}"] += 1
            filePath = '/'.join(m1['relativePath'].split('.')[:-1])
            zuhepath = f"/home/changzexing/RandomMix/{pid}/{pid.lower()}_{vid}_buggy/{filePath}/{linem1}+{linem2}/{flagHOM[m1['relativePath']][f'{linem1}+{linem2}']}"
            # print(zuhepath, filename)
            # break
            if not os.path.exists(zuhepath):
                os.makedirs(zuhepath)
            zuhepath += f"/{filename}"
            # print(zuhepath)
            # path1 = m1["mutFilePath"].replace("/home/fanluxi/pmbfl/mutantsFile", "/home/changzexing/mutantFaultyFileMajor")
            # path2 = m2["mutFilePath"].replace("/home/fanluxi/pmbfl/mutantsFile", "/home/changzexing/mutantFaultyFileMajor")
            path1 = m1["mutFilePath"]
            path2 = m2["mutFilePath"]
            originpath = f"/home/changzexing/d4jclean/{pid}/{vid}b/{m1['relativePath']}"
            # print(vj[0], vj[1], v, vl)
            # print(path1, path2, zuhepath, originpath)
            # break
            command = f"""
java -cp /home/changzexing/d4jscript/codeMerge-1.0-SNAPSHOT.jar:/home/changzexing/d4jscript/javaparser-core-3.24.0.jar:/home/changzexing/d4jscript/java-diff-utils-4.12.jar:/home/changzexing/d4jscript/commons-cli-1.4.jar org.example.Main -pathOri {originpath} -pathM1 {path1} -pathM2 {path2} -pathOutput {zuhepath}
"""
            copyCodeFlag = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            if copyCodeFlag.returncode == 0:
                print(f"{pid}-{vid}-{m1['relativePath']}/{linem1}+{linem2}/{flagHOM[m1['relativePath']][f'{linem1}+{linem2}']}成功")
                print(copyCodeFlag.stdout)
                success += 1
            else:
                print(f"{pid}-{vid}-{m1['relativePath']}/{linem1}+{linem2}/{flagHOM[m1['relativePath']][f'{linem1}+{linem2}']}失败")
                print(copyCodeFlag.stderr)
                fail += 1
        #     break
        # break
        print(f"{pid}-{vid} 完成 成功{success} 失败{fail}")

if __name__ == '__main__':
    init('Chart', 14, 14)
    # init('Time', 1, 27)
    # init('Lang', 1, 65)
    # init('Math', 1, 106)