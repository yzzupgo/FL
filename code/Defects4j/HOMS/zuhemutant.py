import json
import os
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
    successFaultyFileMutant = read_json_file(f"/home/changzexing/successFaultyFileMutant/{pid}.json")
    juzhen = read_json_file(f"/home/changzexing/juzhen/{pid}.json")
    for vid in range(svid, evid + 1):
        print(f"{pid}-{vid}开始")
        success = 0
        fail = 0
        if f"{pid}-{vid}" not in juzhen or f"{pid}-{vid}" not in successFaultyFileMutant:
            continue
        tmpjuzhen = juzhen[f"{pid}-{vid}"]
        tmpsuccess = successFaultyFileMutant[f"{pid}-{vid}"]
        for file in tmpjuzhen:
            filename = file.split('/')[-1]
            for index, value in enumerate(tmpjuzhen[file]):
                for j, vj in enumerate(value):
                    if j < index:
                        continue
                    if str(vj[0]) in tmpsuccess[file] and str(vj[1]) in tmpsuccess[file]:
                      len1 = len(tmpsuccess[file][str(vj[0])])
                      len2 = len(tmpsuccess[file][str(vj[1])])
                      vj1 = tmpsuccess[file][str(vj[0])]
                      vj2 = tmpsuccess[file][str(vj[1])]
                      num = 0
                      if vj[0] == vj[1] and len1 != 1:
                          for i, v in  enumerate(vj1):
                              if num > vj[-1]: break
                              for k, vl in  enumerate(vj2):
                                  if num > vj[-1]: break
                                  if k <= i: continue
                                  num += 1
                                  zuhepath = f"/home/changzexing/zuhemutant/{pid}/{pid.lower()}_{vid}_buggy/{file[:-5]}/{vj[0]}+{vj[1]}/{num}"
                                  if not os.path.exists(zuhepath):
                                    os.makedirs(zuhepath)
                                  zuhepath += f"/{filename}"
                                  path1 = f"/home/changzexing/mutantFaultyFile/{pid}/{pid.lower()}_{vid}_buggy/{file[:-5]}/{vj[0]}/{v}/{filename}"
                                  path2 = f"/home/changzexing/mutantFaultyFile/{pid}/{pid.lower()}_{vid}_buggy/{file[:-5]}/{vj[1]}/{vl}/{filename}"
                                  originpath = f"/home/changzexing/d4jclean/{pid}/{vid}b/{file}"
                                  # print(vj[0], vj[1], v, vl)
                                  # print(path1, path2, zuhepath)
                                  command = f"""
java -cp /home/changzexing/d4jscript/codeMerge-1.0-SNAPSHOT.jar:/home/changzexing/d4jscript/javaparser-core-3.24.0.jar:/home/changzexing/d4jscript/java-diff-utils-4.12.jar:/home/changzexing/d4jscript/commons-cli-1.4.jar org.example.Main -pathOri {originpath} -pathM1 {path1} -pathM2 {path2} -pathOutput {zuhepath}
"""
                                  copyCodeFlag = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                                  if copyCodeFlag.returncode == 0:
                                      print(f"{pid}-{vid}-{filename}-{vj[0]}-{vj[1]}成功")
                                      success += 1
                                  else:
                                      print(f"{pid}-{vid}-{filename}-{vj[0]}-{vj[1]}失败")
                                      fail += 1
                      if vj[0] != vj[1]:
                          for i, v in  enumerate(vj1):
                              if num > vj[-1]: break
                              for k, vl in  enumerate(vj2):
                                  if num > vj[-1]: break
                                  if k < i: continue
                                  num += 1
                                  zuhepath = f"/home/changzexing/zuhemutant/{pid}/{pid.lower()}_{vid}_buggy/{file[:-5]}/{vj[0]}+{vj[1]}/{num}"
                                  if not os.path.exists(zuhepath):
                                    os.makedirs(zuhepath)
                                  zuhepath += f"/{filename}"
                                  path1 = f"/home/changzexing/mutantFaultyFile/{pid}/{pid.lower()}_{vid}_buggy/{file[:-5]}/{vj[0]}/{v}/{filename}"
                                  path2 = f"/home/changzexing/mutantFaultyFile/{pid}/{pid.lower()}_{vid}_buggy/{file[:-5]}/{vj[1]}/{vl}/{filename}"
                                  originpath = f"/home/changzexing/d4jclean/{pid}/{vid}b/{file}"
                                  # print(vj[0], vj[1], v, vl) 
                                  command = f"""
java -cp /home/changzexing/d4jscript/codeMerge-1.0-SNAPSHOT.jar:/home/changzexing/d4jscript/javaparser-core-3.24.0.jar:/home/changzexing/d4jscript/java-diff-utils-4.12.jar:/home/changzexing/d4jscript/commons-cli-1.4.jar org.example.Main -pathOri {originpath} -pathM1 {path1} -pathM2 {path2} -pathOutput {zuhepath}
"""
                                  copyCodeFlag = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                                  if copyCodeFlag.returncode == 0:
                                      print(f"{pid}-{vid}-{filename}-{vj[0]}-{vj[1]}成功")
                                      success += 1
                                  else:
                                      print(f"{pid}-{vid}-{filename}-{vj[0]}-{vj[1]}失败")
                                      fail += 1
                      # print(len1, len2, vj[-1])                   
        # break
        print(f"{pid}-{vid}完成, 成功{success}个,失败{fail}个")

if __name__ == "__main__":
    init("Chart", 7, 26)