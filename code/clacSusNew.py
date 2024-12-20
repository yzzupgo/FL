import ast
import csv
import glob
import json
import os
import subprocess
from tkinter import NO

def read_csv(file_path):
    """
    读取 CSV 文件并返回其中的数据。

    Parameters:
    - file_path (str): CSV 文件的路径。

    Returns:
    - list of dict: CSV 文件中的数据，每行作为一个字典。
    """
    data = []

    try:
        with open(file_path, newline='', encoding='utf-8') as csvfile:
            reader = csv.DictReader(csvfile, skipinitialspace=True)
            for row in reader:
                data.append(row)
        return data
    except FileNotFoundError:
        print(f"Error: File not found at path {file_path}.")
        return None
    except csv.Error as e:
        print(f"Error reading CSV file {file_path}: {e}")
        return None


def find_csv_files(directory):
    """
    查找给定路径下的所有 CSV 文件。

    Parameters:
    - directory (str): 要查找的目录路径。

    Returns:
    - list of str: 匹配的 CSV 文件路径列表。
    """
    csv_files = glob.glob(f"{directory}/*.csv")
    return csv_files


#获取变异体结果json文件
def get_non_txt_files(folder_path):
  files = glob.glob(os.path.join(folder_path, "*"))
  non_txt_files = [os.path.basename(f) for f in files if not f.endswith('.txt') and os.path.isfile(f)]
  return non_txt_files

#找到测试用用例txt及其结果txt
def find_test_txt_files(file_path, txt_name):
  result = subprocess.run(["find", file_path, "-name", txt_name + ".txt"], stdout=subprocess.PIPE)
  txt_files = result.stdout.decode("utf-8").strip().split("\n")
  return txt_files

#逐行读取原始txt
def read_txt_file(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        lines = file.readlines()
    lines = [line.strip() for line in lines]
    return lines

#得到原版本测试用例执行结果：【测试用例：执行结果】
def get_init_test_result(cover_path):
  tests_txt = find_test_txt_files(cover_path, "all_tests")
  tests_result_txt = find_test_txt_files(cover_path, "inVector")
  # print(tests_txt, tests_result_txt)
  test_result = {}
  num = 0
  for item in range(len(tests_txt)):
    test_list = read_txt_file(tests_txt[item]) 
    result_list = read_txt_file(tests_result_txt[item])
    # print(len(test_list), len(result_list), test_list[2], 1 if test_list[2] == 'org.jfree.chart.annotations.junit.CategoryLineAnnotationTests#testCloning' else 0)
    for test_item in range(len(test_list)):
      # print(test_list[test_item], result_list[test_item])
      test_result[str(test_list[test_item])] = result_list[test_item]
      if str(result_list[test_item]) =='1':
         num += 1
      # if test_item <= 10:
      #   print(test_list[test_item], result_list[test_item])
  print("1的数量为{}".format(num))
  return test_result

#得到变异体测试用例
def mutant_test(filename):
    with open(filename) as f:
        lines = f.readlines()

    processed_lines = []
    for line in lines:
        line = line.strip()
        left_idx = line.find("(")
        right_idx = line.find(")")
        if left_idx != -1 and right_idx != -1:
            class_name = line[:left_idx].split(".")[-1]
            method_name = line[left_idx+1:right_idx].replace(",", "#")
            processed_line = f"{method_name}#{class_name}"
            processed_lines.append(processed_line)

    return processed_lines


def read_json_file(file_path):
    """
    读取JSON文件并返回解析后的数据

    Parameters:
    - file_path (str): JSON文件的路径

    Returns:
    - data (dict): 解析后的JSON数据
    """
    try:
        with open(file_path, 'r') as file:
            data = json.load(file)
            return data
    except FileNotFoundError:
        print(f"文件 '{file_path}' 未找到。")
        return None
    except json.JSONDecodeError as e:
        print(f"解析JSON时发生错误: {e}")
    except Exception as e:
        print(f"发生错误: {e}")

def Ochiai(kf, failed, killed):
    if failed == 0 or killed == 0:
      return 0
    # if kf == 0:
    #    kf = 1
    return kf / ((failed * killed) ** 0.5)

def Dstar(kf, kp, nf):
  if kp + nf == 0:
    return 0
  return (kf * kf) / (kp + nf)

def gp13(Akf, Anf, Akp, Anp):
    if (2 * Akp + Akf) == 0:
        return 0
    return Akf + (Akf / (2 * Akp + Akf))

def Jaccard(kf, kp, nf):
   if (kf + nf + kp == 0):
      return 0
   return kf / (kf + nf + kp)

def Tarantula(kf, kp, nf, np):
   if kf + kp == 0 or kf + nf == 0:
      return 0
   if kf + nf != 0 and kp + np == 0:
      return 1
   return (kf / (kf + nf)) / (kf / (kf + nf) + kp / (kp + np))

def Op2(kf, kp, np):
   return kf - kp / (kp + np + 1)

def init(init_path, mutant_path, mutant_result_path, vid, pid, susName):
  print(init_path, os.path.exists(init_path))
  print(mutant_path, os.path.exists(mutant_path))
  print(mutant_path + '/all_tests', os.path.exists(mutant_path + '/all_tests'))
  # if not os.path.exists(init_path) or not os.path.exists(mutant_path) or not os.path.exists(mutant_path + '/all_tests'):
  if not os.path.exists(init_path):
    print("路径缺失")
    return 
  print(pid, vid)
  initresult = get_init_test_result(init_path) #Gzotar工具获取的所有测试用例
  # mutant_result = mutant_test(mutant_path + '/all_tests') #执行源程序test后获得的测试用例
  # print(result) 
  # print(len(initresult), len(mutant_result))
  # mixtest = [] #给 initresult 和 mutant_result求一个交集
  mixtest = initresult
  # 记录major变异体和处理成json之后的映射关系
  mapmajor = {}
  # 处理major的变异算子类型开始
  # mutantInfo = {}
  # # linemap 记录major变异体的行号，某个版本存在多个文件时由于会存在相同的行号，行号会顺延，例如A文件100行有4个，那么B文件的100行变异体从5开始
  # linemap = {}
  # # item_path_major 是major的变异体信息路径
  # # item_path_major = f"/home/changzexing/mutant_result_faulty_file_major/{pid}/{vid}b/muInfo.json"
  # item_path_major = f"/home/changzexing/mutant_result_faulty_file_major/{pid}/{vid}b/mutantInfo.json"
  # majormutantresultinfo = read_json_file(item_path_major)
  # if majormutantresultinfo is None:
  #    return
  # for item_major in majormutantresultinfo:
  #     num = item_major['linenum']
  #     # print(type(num))
  #     if num in linemap:
  #        linemap[num] += 1
  #     else:
  #        linemap[num] = 1
  #     # 原本的格式，只有Math 和 Closure的1-33能用
  #     # fileName = f"{item_major['relativePath'].split('.')[0].replace('/', '-')}-{num}"
  #     # 租的服务器的格式
  #     fileName = '-'.join(item_major['relativePath'].split('.')[0].split('/')[6:]) + f'-{num}'
  #     if fileName not in mutantInfo:
  #       mutantInfo[fileName] = {}
  #     mutantInfo[fileName][str(linemap[num])] = {
  #       'mut_operator': item_major["typeOp"],
  #       'bianhao': linemap[num]
  #     }
  #     mapmajor[f"{fileName}-{linemap[num]}"] = item_major
  # print(linemap)
  # 处理major的变异算子类型结束
  killJuzhen = {}
  killJuzhenMaxSus = {}
  susdic_Oc = {}
  susdic_Ds = {}
  susdic_Ja = {}
  susdic_Op = {}
  susdic_Ta = {}
  susdic_Gp = {}
  maxMbertMutantInfo_Oc = {}
  maxMbertMutantInfo_Ds = {}
  maxMbertMutantInfo_Ja = {}
  maxMbertMutantInfo_Op = {}
  maxMbertMutantInfo_Ta = {}
  maxMbertMutantInfo_Gp = {}
  sum_susdic_Oc = {}
  sum_susdic_Ds = {}
  sum_susdic_Ja = {}
  sum_susdic_Op = {}
  sum_susdic_Ta = {}
  sum_susdic_Gp = {}
  num_susdic_Oc = {}
  num_susdic_Ds = {}
  num_susdic_Ja = {}
  num_susdic_Op = {}
  num_susdic_Ta = {}
  num_susdic_Gp = {}
  set_susdic_Oc = {}
  set_susdic_Ds = {}
  set_susdic_Ja = {}
  set_susdic_Op = {}
  set_susdic_Ta = {}
  set_susdic_Gp = {}
  # 取交集
  # for key in initresult.keys():
  #    if key in mutant_result:
  #       mixtest.append(key)
  # for item in mutant_result:
  #    if item in initresult and item not in mixtest:
  #       mixtest.append(key)
  # print(len(result), len(mutant_result), len(mixtest))
  mutant_result_json = get_non_txt_files(mutant_result_path) # 变异体执行结果
  # print(len(mutant_result_json))
  for item in mutant_result_json: 
    line_name = "-".join(item.split("-")[:-1])
    # print(line_name, item)
    # susdic_Oc代表每行的怀疑度，sum_susdic_Oc表示该行怀疑度的总和，num_susdic表示该行变异体的总数
    susdic_Oc[line_name] = -1
    susdic_Ds[line_name] = -1
    susdic_Ja[line_name] = -1
    susdic_Op[line_name] = -1
    susdic_Ta[line_name] = -1
    susdic_Gp[line_name] = -1
    sum_susdic_Oc[line_name] = 0
    sum_susdic_Ds[line_name] = 0
    sum_susdic_Ja[line_name] = 0
    sum_susdic_Op[line_name] = 0
    sum_susdic_Ta[line_name] = 0
    sum_susdic_Gp[line_name] = 0
    num_susdic_Oc[line_name] = 0
    num_susdic_Ds[line_name] = 0
    num_susdic_Ja[line_name] = 0
    num_susdic_Op[line_name] = 0
    num_susdic_Ta[line_name] = 0
    num_susdic_Gp[line_name] = 0

  # 遍历每个变异体，获取每个变异体的怀疑度
  for item in mutant_result_json:
    print(len(mutant_result_json))
    juzhenkey = item.split('.')[0]
    killJuzhen[juzhenkey] = [] # 高阶所需要的文件
    # itme=>src-main-java-org-apache-commons-math3-fraction-Fraction-188-17.json
    print(line_name)
    line_name = "-".join(item.split("-")[:-1]) # src-main-java-org-apache-commons-math3-fraction-Fraction-188 188是这个文件对应的行号
    # item_path 是μbert的变异体路径
    item_path = f"/home/changzexing/mutantFaultyFile/{pid}/{pid.lower()}_{vid}_buggy/{'/'.join(item[:-5].split('-')[:-1])}"
    mutantNo = item[:-5].split('-')[-1]
    # 处理μbert的变异算子类型开始
    # mutantInfo = {} # 存储当前行所有变异体的变异算子信息 
    # mutantInfo[line_name] = {}
    # csvfile = find_csv_files(item_path)
    # for csvfile_item in csvfile: 
    #   csvdata = read_csv(csvfile_item)
    #   # print('csv', csvdata)
    #   for csvdata_item in csvdata:
    #      mutantInfo[line_name][csvdata_item['id']] = {
    #         'mut_operator': csvdata_item["mut_operator"],
    #         'orig_token': csvdata_item['orig_token'],
    #         'pred_token': csvdata_item['pred_token'],
    #         'bainhao': csvdata_item['id']
    #      }
    # 处理μbert的变异算子类型结束
    # print(csvfile, type(csvfile))
    # print(item_path, item)
    # print(mutantInfo) 
    with open(mutant_result_path + '/' + item) as f:
      data = json.load(f)
      kf = nf = kp = np = 0
      # 遍历每个测试用例，统计参数杀死信息
      for test_item in mixtest:
        # pass
        if str(initresult[test_item]) == '0':
          # if test_item.split('#')[0] in data:
          #    kp +=1
          # elif test_item in data:
          if test_item in data:
             kp += 1
             killJuzhen[juzhenkey].append(1)
          else:
             np += 1
             killJuzhen[juzhenkey].append(0)
        # fail
        elif str(initresult[test_item]) == '1':
          # print(vid + 'b', item, test_item, initresult[test_item])
          with open('/home/changzexing/failingTestOutput/' + pid + '/'+ vid + 'b/failing_tests.json') as f:
            faildata = json.load(f) # 数据说明，faildata是原始文件的type1、2、3、4错误信息，会有若干个错误的测试类，组成一个大的json对象，test_item是每个测试类的名称
            # print(data)
            # input()
            # 这里如果测试的项目不在原始程序的错误测试类中,直接跳过
            # if test_item not in faildata:
            #    continue
            # if test_item.split('#')[0] in data:
            #    nf += 1
              # if data[test_item.split('#')[0]].split('at')[0] == faildata[test_item].split('at')[0]:
              # if data[test_item.split('#')[0]] == faildata[test_item]:  
              # if data[test_item.split('#')[0]]['type3'] == faildata[test_item].split('at')[0]: 
              # if data[test_item.split('#')[0]]['type4'] == faildata[test_item]: 
              #   nf += 1
              # else:
              #   kf += 1
            if test_item in data:
              #  nf += 1
              #  killJuzhen[juzhenkey].append(0)
              # 下面全注释表示type1
              # if data[test_item].split('at')[0] == faildata[test_item].split('at')[0]: 
              # if data[test_item] == faildata[test_item]:pytho
              # if data[test_item]['type3'] == faildata[test_item]['type3']:
              # print(test_item, data[test_item]['type3'], faildata[test_item])
              if test_item in faildata and data[test_item]['type3'] == faildata[test_item]['type3']:
                nf += 1
                killJuzhen[juzhenkey].append(0)
              else:
                kf += 1
                killJuzhen[juzhenkey].append(1)
            else:
              kf += 1
              killJuzhen[juzhenkey].append(1)
  
      # 得到中间文件：
      # myfile = '/home/changzexing/susMaxMajorType1/' + pid + '/' + vid + 'b'
      # if not os.path.exists(myfile):
      #   os.makedirs(myfile) 
      # with open(myfile + '/b.txt', mode='a', newline='') as file:
      #   print(item, kf, kp, nf, np, file=file)
      Oc = Ochiai(kf, kf + nf, kf + kp)
      Ds = Dstar(kf, kp, nf)
      Ja = Jaccard(kf, kp, nf)
      Op = Op2(kf, kp, np)
      Ta = Tarantula(kf, kp, nf, np)
      Gp = gp13(kf, nf, kp, np)

      set_susdic_Oc[juzhenkey] = Oc
      set_susdic_Ds[juzhenkey] = Ds
      set_susdic_Ja[juzhenkey] = Ja
      set_susdic_Op[juzhenkey] = Op
      set_susdic_Ta[juzhenkey] = Ta
      set_susdic_Gp[juzhenkey] = Gp
      # avg==其实是sum
      sum_susdic_Oc[line_name] = sum_susdic_Oc[line_name] + Oc
      sum_susdic_Ds[line_name] = sum_susdic_Ds[line_name] + Ds
      sum_susdic_Ja[line_name] = sum_susdic_Ja[line_name] + Ja
      sum_susdic_Op[line_name] = sum_susdic_Op[line_name] + Op
      sum_susdic_Ta[line_name] = sum_susdic_Ta[line_name] + Ta
      sum_susdic_Gp[line_name] = sum_susdic_Gp[line_name] + Gp
      num_susdic_Oc[line_name] = num_susdic_Oc[line_name] + 1
      num_susdic_Ds[line_name] = num_susdic_Ds[line_name] + 1
      num_susdic_Ja[line_name] = num_susdic_Ja[line_name] + 1
      num_susdic_Op[line_name] = num_susdic_Op[line_name] + 1
      num_susdic_Ta[line_name] = num_susdic_Ta[line_name] + 1
      num_susdic_Gp[line_name] = num_susdic_Gp[line_name] + 1
      # if Oc > susdic_Oc[line_name]:
      #   #  print(line_name, mutantInfo[line_name], mutantNo)
      #    maxMbertMutantInfo_Oc[line_name] = mutantInfo[line_name][mutantNo]
      # if Ds > susdic_Ds[line_name]:
      #    maxMbertMutantInfo_Ds[line_name] = mutantInfo[line_name][mutantNo]
      # if Ja > susdic_Ja[line_name]:
      #    maxMbertMutantInfo_Ja[line_name] = mutantInfo[line_name][mutantNo]
      # if Op > susdic_Op[line_name]:
      #    maxMbertMutantInfo_Op[line_name] = mutantInfo[line_name][mutantNo]
      # if Ta > susdic_Ta[line_name]:
      #    maxMbertMutantInfo_Ta[line_name] = mutantInfo[line_name][mutantNo]
      # if Gp > susdic_Gp[line_name]:
      #    maxMbertMutantInfo_Gp[line_name] = mutantInfo[line_name][mutantNo]
      susdic_Oc[line_name] = max(susdic_Oc[line_name], Ochiai(kf, kf + nf, kf + kp))
      susdic_Ds[line_name] = max(susdic_Ds[line_name], Dstar(kf, kp, nf))
      susdic_Ja[line_name] = max(susdic_Ja[line_name], Jaccard(kf, kp, nf))
      susdic_Op[line_name] = max(susdic_Op[line_name], Op2(kf, kp, np))
      susdic_Ta[line_name] = max(susdic_Ta[line_name], Tarantula(kf, kp, nf, np))
      susdic_Gp[line_name] = max(susdic_Gp[line_name], gp13(kf, nf, kp, np))
      # susdic_Oc[line_name].append(Ochiai(kf, kf + nf, kf + kp))
      # susdic_Ds[line_name].append(Dstar(kf, kp, nf))
      # susdic_Ja[line_name].append(Jaccard(kf, kp, nf))
    #   print('{} {}b {} 计算完成'.format(pid, vid, item))
      # print(maxMbertMutantInfo_Oc)
      # break
  # 求average作为muse的结果：
  if 'Muse' in susName:
    for tmp in susdic_Oc.keys():
      susdic_Oc[tmp] = sum_susdic_Oc[tmp] / num_susdic_Oc[tmp]
    for tmp in susdic_Ds.keys():
      susdic_Ds[tmp] = sum_susdic_Ds[tmp] / num_susdic_Ds[tmp]
    for tmp in susdic_Ja.keys():
      susdic_Ja[tmp] = sum_susdic_Ja[tmp] / num_susdic_Ja[tmp]
    for tmp in susdic_Op.keys():
      susdic_Op[tmp] = sum_susdic_Op[tmp] / num_susdic_Op[tmp]
    for tmp in susdic_Ta.keys():
      susdic_Ta[tmp] = sum_susdic_Ta[tmp] / num_susdic_Ta[tmp]
    for tmp in susdic_Gp.keys():
      susdic_Gp[tmp] = sum_susdic_Gp[tmp] / num_susdic_Gp[tmp]    
  
  
  # 求max-average
  # for tmp in susdic_Oc.keys():
  #   susdic_Oc[tmp] = susdic_Oc[tmp] - sum_susdic_Oc[tmp] / num_susdic_Oc[tmp]
  # for tmp in susdic_Ds.keys():
  #   susdic_Ds[tmp] = susdic_Ds[tmp] - sum_susdic_Ds[tmp] / num_susdic_Ds[tmp]
  # for tmp in susdic_Ja.keys():
  #   susdic_Ja[tmp] = susdic_Ja[tmp] - sum_susdic_Ja[tmp] / num_susdic_Ja[tmp]
  # for tmp in susdic_Op.keys():
  #   susdic_Op[tmp] = susdic_Op[tmp] - sum_susdic_Op[tmp] / num_susdic_Op[tmp]
  # for tmp in susdic_Ta.keys():
  #   susdic_Ta[tmp] = susdic_Ta[tmp] - sum_susdic_Ta[tmp] / num_susdic_Ta[tmp]
  # for tmp in susdic_Gp.keys():
  #   susdic_Gp[tmp] = susdic_Gp[tmp] - sum_susdic_Gp[tmp] / num_susdic_Gp[tmp]

  # 将怀疑度列表排序
  sorted_dict_Oc = dict(sorted(susdic_Oc.items(), key=lambda item: item[1], reverse=True))
  sorted_dict_Ds = dict(sorted(susdic_Ds.items(), key=lambda item: item[1], reverse=True))
  sorted_dict_Ja = dict(sorted(susdic_Ja.items(), key=lambda item: item[1], reverse=True))
  sorted_dict_Op = dict(sorted(susdic_Op.items(), key=lambda item: item[1], reverse=True))
  sorted_dict_Ta = dict(sorted(susdic_Ta.items(), key=lambda item: item[1], reverse=True))
  sorted_dict_Gp = dict(sorted(susdic_Gp.items(), key=lambda item: item[1], reverse=True))

  print(sorted_dict_Ds)
  
  # myfile = '/home/changzexing/susMaxMbertMutantInfo/' + pid + '/' + vid + 'b'
  # myfile = '/home/changzexing/susMaxMajorMutantInfo/' + pid + '/' + vid + 'b'
  # myfile = '/home/changzexing/MaxMajorKillJuzhen/' + pid + '/' + vid + 'b'
  # myfile = '/home/changzexing/MaxMbertKillJuzhen/' + pid + '/' + vid + 'b'
  # 将排序后的字典写入csv文件,选择对应的文件
  # myfile = '/home/changzexing/susMaxMbertHebingType1/' + pid + '/' + vid + 'b'
  # myfile = '/home/changzexing/susMaxSubAvgMbertType1/' + pid + '/' + vid + 'b'
  # myfile = '/home/changzexing/susMaxSubAvgMbertType3/' + pid + '/' + vid + 'b'
  # myfile = '/home/changzexing/susMaxSubAvgMbertType4/' + pid + '/' + vid + 'b'
  # myfile = '/home/changzexing/susMaxSubAvgMajorType1/' + pid + '/' + vid + 'b'
  # myfile = '/home/changzexing/susMaxSubAvgMajorType3/' + pid + '/' + vid + 'b'
  # myfile = '/home/changzexing/susMaxSubAvgMajorType4/' + pid + '/' + vid + 'b'
  # myfile = '/home/changzexing/susMaxMbertType1/' + pid + '/' + vid + 'b'
  # myfile = '/home/changzexing/susMaxMbertType3/' + pid + '/' + vid + 'b'
  # myfile = '/home/changzexing/susMaxMbertType4/' + pid + '/' + vid + 'b'
  # myfile = '/home/changzexing/susMaxMajorType1/' + pid + '/' + vid + 'b'
  # myfile = '/home/changzexing/susMaxMajorType3/' + pid + '/' + vid + 'b'
  # myfile = '/home/changzexing/susMaxMajorType4/' + pid + '/' + vid + 'b'
  # myfile = '/home/changzexing/susMaxSFCluShanjian/' + pid + '/' + vid + 'b'
  # myfile = '/home/changzexing/susMaxSFClu/' + pid + '/' + vid + 'b'
  mp = {
      'susMaxMajorType3': 'Major',
      'susMaxMbertType3': 'Mbert'
  }
  # 获取杀死矩阵
  # killJuzhenFile = f'/home/changzexing/Max{mp[susName]}KillJuzhen/{pid}/{vid}b'
  # killJuzhenFile = f'/home/changzexing/MaxMbertKillJuzhen/{pid}/{vid}b'
  # if not os.path.exists(killJuzhenFile):
  #   os.makedirs(killJuzhenFile)
  # with open(f"{killJuzhenFile}/killJuzhen.json", 'w') as f:
  #   json.dump(killJuzhen, f)
  # return 
  # 存储major原始结果和处理成json之后的映射关系
  # with open(f"/home/changzexing/MaxMajorKillJuzhen/{pid}/{vid}b/mapmajor.json", 'w') as f:
  #   json.dump(mapmajor, f)
  # return

  # print(set_susdic_Oc)

  # 获取每个变异体的怀疑度
  # AllMutantSusFile = f'/home/changzexing/All{mp[susName]}MutantSus/{pid}/{vid}b'
  # if not os.path.exists(AllMutantSusFile):
  #   os.makedirs(AllMutantSusFile)
  # with open(f"{AllMutantSusFile}/Ochiai.json", 'w') as f:
  #   json.dump(set_susdic_Oc, f)
  # with open(f"{AllMutantSusFile}/Dstar.json", 'w') as f:
  #   json.dump(set_susdic_Ds, f)
  # with open(f"{AllMutantSusFile}/Jaccard.json", 'w') as f:
  #   json.dump(set_susdic_Ja, f)
  # with open(f"{AllMutantSusFile}/Op2.json", 'w') as f:
  #   json.dump(set_susdic_Op, f)
  # with open(f"{AllMutantSusFile}/Tarantula.json", 'w') as f:
  #   json.dump(set_susdic_Ta, f)
  # with open(f"{AllMutantSusFile}/Gp13.json", 'w') as f:
  #   json.dump(set_susdic_Gp, f)
  # return
    
  # 处理变异算子
  # mutantInfoFile = f'/home/changzexing/susMaxMajorMutantInfo/{pid}/{vid}b'
  # # mutantInfoFile = f'/home/changzexing/susMaxMbertMutantInfo/{pid}/{vid}b'
  # if not os.path.exists(mutantInfoFile):
  #   os.makedirs(mutantInfoFile)
  # with open(f"{mutantInfoFile}/Ochiai.json", 'w') as f:
  #   json.dump(maxMbertMutantInfo_Oc, f)
  # with open(f"{mutantInfoFile}/Dstar.json", 'w') as f:
  #   json.dump(maxMbertMutantInfo_Ds, f)
  # with open(f"{mutantInfoFile}/Jaccard.json", 'w') as f:
  #   json.dump(maxMbertMutantInfo_Ja, f)
  # with open(f"{mutantInfoFile}/Op2.json", 'w') as f:
  #   json.dump(maxMbertMutantInfo_Op, f)
  # with open(f"{mutantInfoFile}/Tarantula.json", 'w') as f:
  #   json.dump(maxMbertMutantInfo_Ta, f)
  # with open(f"{mutantInfoFile}/Gp13.json", 'w') as f:
  #   json.dump(maxMbertMutantInfo_Gp, f)
  
  # return
  # 获取怀疑度
  myfile = '/home/changzexing/' + susName + '/' + pid + '/' + vid + 'b'
  # print(myfile)
  if not os.path.exists(myfile):
    os.makedirs(myfile)
  with open(myfile + '/Ochiai.csv', mode='w', newline='') as file:
      writer = csv.writer(file)
      for key, value in sorted_dict_Oc.items():
          writer.writerow([key, value])
  with open(myfile + '/Dstar.csv', mode='w', newline='') as file:
      writer = csv.writer(file)
      for key, value in sorted_dict_Ds.items():
          writer.writerow([key, value])
  with open(myfile + '/Jaccard.csv', mode='w', newline='') as file:
      writer = csv.writer(file)
      for key, value in sorted_dict_Ja.items():
          writer.writerow([key, value])
  with open(myfile + '/Op2.csv', mode='w', newline='') as file:
      writer = csv.writer(file)
      for key, value in sorted_dict_Op.items():
          writer.writerow([key, value])
  with open(myfile + '/Tarantula.csv', mode='w', newline='') as file:
      writer = csv.writer(file)
      for key, value in sorted_dict_Ta.items():
          writer.writerow([key, value])
  with open(myfile + '/Gp13.csv', mode='w', newline='') as file:
      writer = csv.writer(file)
      for key, value in sorted_dict_Gp.items():
          writer.writerow([key, value])
     
def start(pid, svid, evid, name, susName):
   for i in range(svid, evid + 1):
      # if str(i) == '2':
      #     continue
      # list_ignore = [5, 6, 13, 14, 15, 16, 17, 19, 36, 44, 54, 59, 71, 73, 74, 99]
      # if i in list_ignore:
      #    continue
      init(f'/home/changzexing/d4jbasecover/{pid}/{i}b', f'/home/changzexing/d4jclean/{pid}/{i}b', f'/home/changzexing/{name}/{pid}/{i}b', str(i), pid, susName)
      
      # init('/home/changzexing/d4jbasecover/' + pid + '/' + str(i) + 'b', '/home/changzexing/d4jclean/' + pid + '/' + str(i) + 'b', '/home/changzexing/mutant_result_faulty_file_major_json/' + pid + '/' + str(i) + 'b', str(i), pid)


def chushihua(name, susName):
    start('Chart', 1, 26, name, susName)
    # start('Time', 1, 27, name, susName)
    # start('Lang', 1, 65, name, susName)
    # start('Math', 1, 106, name, susName)
    # start('Mockito', 1, 38, name, susName)
    # start('Closure', 1, 133, name, susName)
    # start('Cli', 1, 39, name, susName)
    # start('Codec', 1, 18, name, susName)
    # start('Csv', 1, 16, name, susName)
    # start('Gson', 1, 18, name, susName)
    # start('JacksonXml', 1, 6, name, susName)

if __name__ == '__main__':
    # 1. 确定309/321文件位置
    # 2. 统一183和257的类型，要求一样，例如183是type3.257也需要是type3
    # 3. 如果聚合方式为max-avg，则需要打开229-240行，注释时是max的情况
    # chushihua('mutant_result_faulty_file_SFClu_shanjian_json')
  #  chushihua('mutant_result_faulty_file_SFClu_json', 'susMaxSFCluQuanjiHOM')
  #  chushihua('mutant_result_faulty_file_Major_Last2First_json', 'susMaxMajorLast2FirstHOM')
  #  chushihua('mutant_result_faulty_file_Major_RandomMix_json', 'susMaxMajorRandomMixHOM')
  #  chushihua('mutant_result_faulty_file_Last2First_MajorMbert_json', 'susMaxMajorMbertLast2FirstHOM')
  #  chushihua('mutant_result_faulty_file_RandomMix_MajorMbert_json', 'susMaxMajorMbertRandomMixHOM')
  #  chushihua('mutant_result_faulty_file_Major_Mbert_SFClu_shanjian_json', 'susMaxMajorMbertSFCluHOM')
  # chushihua('mutant_result_faulty_file_major_json', 'MusesusMaxMajorType3')
  # chushihua('mutant_result_faulty_file_json', 'MusesusMaxMbertType3')
  # chushihua('mutant_result_faulty_file_Major_SFClu_json', 'MusesusMaxMajorSFCluHOMType3')
  # chushihua('mutant_result_faulty_file_Major_Mbert_SFClu_json', 'MusesusMaxMajorMbertSFCluHOMType3')
  # chushihua('mutant_result_faulty_file_Major_Last2First_json', 'MusesusMaxMajorLast2FirstType3')
  # chushihua('mutant_result_faulty_file_Major_RandomMix_json', 'MusesusMaxMajorRandomMixType3')
  # chushihua('mutant_result_faulty_file_Major_SamelineHOM_json', 'MusesusMaxMajorSamelineHOMType3')
  # chushihua('mutant_result_faulty_file_Major_Mbert_SamelineHOM_json', 'MusesusMaxMajorMbertSamelineHOMType3')
  
  # chushihua('mutant_result_faulty_file_major_json', 'susMaxMajorType3')
  # chushihua('mutant_result_faulty_file_json', 'susMaxMbertType3')
  # chushihua('mutant_result_faulty_file_Major_SFClu_json', 'susMaxMajorSFCluHOMType3')
  # chushihua('mutant_result_faulty_file_Major_Mbert_SFClu_json', 'susMaxMajorMbertSFCluHOMType3')
  # chushihua('mutant_result_faulty_file_Major_Last2First_json', 'susMaxMajorLast2FirstType3')
  # chushihua('mutant_result_faulty_file_Major_RandomMix_json', 'susMaxMajorRandomMixType3')
  # chushihua('mutant_result_faulty_file_Major_SamelineHOM_json', 'susMaxMajorSamelineHOMType3')
  # chushihua('mutant_result_faulty_file_Major_Mbert_SamelineHOM_json', 'susMaxMajorMbertSamelineHOMType3')
  #  chushihua('mutant_result_faulty_file_json', 'susMaxMbertType3')
  #  chushihua('mutant_result_faulty_file_Major_Mbert_SFClu_json', 'susMaxMajorMbertSFCluHOMType3')
  #  chushihua('mutant_result_faulty_file_Major_SFClu_FaultyLine_json', 'susMaxSFCluShanjianFaultyLineHOM')
  #  chushihua('mutant_result_faulty_file_Major_SFClu_json', 'susMaxMajorSFCluHOM')

  chushihua('mutant_result_faulty_file_Major_zuhemutantSFDis_json', 'susMaxMajorSFDisHOMType3')
  #  chushihua('mutant_result_faulty_file_SFClu_json')
    # start('Chart', 1, 26, '')
    # start('Time', 1, 27)
    # start('Lang', 1, 65)
    # start('Math', 1, 106)
    # pid = 'Time'
    # svid = 1
    # evid = 27
    # pid = 'Chart'
    # svid = 1
    # evid = 26
    # pid = 'Math'
    # svid = 1
    # evid = 106
    # pid = 'Lang'
    # svid = 1
    # evid = 65
    # for i in range(svid, evid + 1):
      # if str(i) == '2':
      #     continue
      # list_ignore = [5, 6, 13, 14, 15, 16, 17, 19, 36, 44, 54, 59, 71, 73, 74, 99]
      # if i in list_ignore:
      #    continue
      # init('/home/changzexing/d4jbasecover/' + pid + '/' + str(i) + 'b', '/home/changzexing/d4jclean/' + pid + '/' + str(i) + 'b', '/home/changzexing/mutant_result_faulty_file_json/' + pid + '/' + str(i) + 'b', str(i), pid)
      
      # init('/home/changzexing/d4jbasecover/' + pid + '/' + str(i) + 'b', '/home/changzexing/d4jclean/' + pid + '/' + str(i) + 'b', '/home/changzexing/mutant_result_faulty_file_major_json/' + pid + '/' + str(i) + 'b', str(i), pid)
      # break

    # pid = 'Time'
    # svid = 1
    # evid = 27

    # for i in range(svid, evid + 1):
    #   # if str(i) == '2':
    #   #     continue
    #   # 1.
    #   init('/home/changzexing/d4jbasecover/' + pid + '/' + str(i) + 'b', '/home/changzexing/d4jclean/' + pid + '/' + str(i) + 'b', '/home/changzexing/mutant_result_faulty_file_json/' + pid + '/' + str(i) + 'b', str(i), pid)
    #   # init('/home/changzexing/d4jbasecover/' + pid + '/' + str(i) + 'b', '/home/changzexing/d4jclean/' + pid + '/' + str(i) + 'b', '/home/changzexing/mutant_result_faulty_file_major_json/' + pid + '/' + str(i) + 'b', str(i), pid)
    #   # break

