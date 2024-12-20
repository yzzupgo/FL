import csv
import glob
import json
import os
import subprocess

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

def init(init_path, mutant_path, mutant_result_path, vid, pid):
  if not os.path.exists(init_path) or not os.path.exists(mutant_path) or not os.path.exists(mutant_path + '/all_tests'):
    return 
  print(pid, vid)
  initresult = get_init_test_result(init_path) #原始测试用例执行结果
  mutant_result = mutant_test(mutant_path + '/all_tests') #变异体测试用例
  # print(result)
  # print(mutant_result)
  # 
  mixtest = []
  susdic_Oc = {}
  susdic_Ds = {}
  susdic_Ja = {}
  susdic_Op = {}
  susdic_Ta = {}
  susdic_Gp = {}
  avg_susdic_Oc = {}
  avg_susdic_Ds = {}
  avg_susdic_Ja = {}
  avg_susdic_Op = {}
  avg_susdic_Ta = {}
  avg_susdic_Gp = {}
  num_susdic_Oc = {}
  num_susdic_Ds = {}
  num_susdic_Ja = {}
  num_susdic_Op = {}
  num_susdic_Ta = {}
  num_susdic_Gp = {}
  for key in initresult.keys():
     if key in mutant_result:
        mixtest.append(key)
  for item in mutant_result:
     if item in initresult and item not in mixtest:
        mixtest.append(key)
  # print(len(result), len(mutant_result), len(mixtest))
  mutant_result_json = get_non_txt_files(mutant_result_path) # 变异体执行结果
  for item in mutant_result_json:
    line_name = "-".join(item.split("-")[:-1])
    susdic_Oc[line_name] = 0
    susdic_Ds[line_name] = 0
    susdic_Ja[line_name] = 0
    susdic_Op[line_name] = 0
    susdic_Ta[line_name] = 0
    susdic_Gp[line_name] = 0
    avg_susdic_Oc[line_name] = 0
    avg_susdic_Ds[line_name] = 0
    avg_susdic_Ja[line_name] = 0
    avg_susdic_Op[line_name] = 0
    avg_susdic_Ta[line_name] = 0
    avg_susdic_Gp[line_name] = 0
    num_susdic_Oc[line_name] = 0
    num_susdic_Ds[line_name] = 0
    num_susdic_Ja[line_name] = 0
    num_susdic_Op[line_name] = 0
    num_susdic_Ta[line_name] = 0
    num_susdic_Gp[line_name] = 0

  # 遍历变异体执行结果的每个测试类,获取每个变异体的怀疑度
  for item in mutant_result_json:
    # itme=>src-main-java-org-apache-commons-math3-fraction-Fraction-188-17.json
    line_name = "-".join(item.split("-")[:-1]) # src-main-java-org-apache-commons-math3-fraction-Fraction-188 188是这个文件对应的行号
    with open(mutant_result_path + '/' + item) as f:
      data = json.load(f)
      kf = nf = kp = np = 0
      for test_item in mixtest:
        # pass
        if str(initresult[test_item]) == '0':
          # if test_item.split('#')[0] in data:
          #    kp +=1
          # elif test_item in data:
          if test_item in data:
             kp += 1
          else:
             np += 1
        # fail
        elif str(initresult[test_item]) == '1':
          # print(vid + 'b', item, test_item, initresult[test_item])
          with open('/home/changzexing/failingTestOutput/' + pid + '/'+ vid + 'b/failing_tests.json') as f:
            faildata = json.load(f) # 数据说明，faildata是原始文件的type1、2、3、4错误信息，会有若干个错误的测试类，组成一个大的json对象，test_item是每个测试类的名称
            # print(data)
            # input()
            # 这里如果测试的项目不在原始程序的错误测试类中,直接跳过
            if test_item not in faildata:
               continue
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
              # 下面全注释表示type1
              # if data[test_item].split('at')[0] == faildata[test_item].split('at')[0]: 
              # if data[test_item] == faildata[test_item]:pytho
              # if data[test_item]['type3'] == faildata[test_item]['type3']:
              if data[test_item]['type3'] == faildata[test_item]['type3']:
                nf += 1
              else:
                kf += 1
            else:
              kf += 1
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
      # avg==其实是sum
      avg_susdic_Oc[line_name] = avg_susdic_Oc[line_name] + Oc
      avg_susdic_Ds[line_name] = avg_susdic_Ds[line_name] + Ds
      avg_susdic_Ja[line_name] = avg_susdic_Ja[line_name] + Ja
      avg_susdic_Op[line_name] = avg_susdic_Op[line_name] + Op
      avg_susdic_Ta[line_name] = avg_susdic_Ta[line_name] + Ta
      avg_susdic_Gp[line_name] = avg_susdic_Gp[line_name] + Gp
      num_susdic_Oc[line_name] = num_susdic_Oc[line_name] + 1
      num_susdic_Ds[line_name] = num_susdic_Ds[line_name] + 1
      num_susdic_Ja[line_name] = num_susdic_Ja[line_name] + 1
      num_susdic_Op[line_name] = num_susdic_Op[line_name] + 1
      num_susdic_Ta[line_name] = num_susdic_Ta[line_name] + 1
      num_susdic_Gp[line_name] = num_susdic_Gp[line_name] + 1
      susdic_Oc[line_name] = max(susdic_Oc[line_name], Ochiai(kf, kf + nf, kf + kp))
      susdic_Ds[line_name] = max(susdic_Ds[line_name], Dstar(kf, kp, nf))
      susdic_Ja[line_name] = max(susdic_Ja[line_name], Jaccard(kf, kp, nf))
      susdic_Op[line_name] = max(susdic_Op[line_name], Op2(kf, kp, np))
      susdic_Ta[line_name] = max(susdic_Ta[line_name], Tarantula(kf, kp, nf, np))
      susdic_Gp[line_name] = max(susdic_Gp[line_name], gp13(kf, nf, kp, np))
      # susdic_Oc[line_name].append(Ochiai(kf, kf + nf, kf + kp)) 
      # susdic_Ds[line_name].append(Dstar(kf, kp, nf))
      # susdic_Ja[line_name].append(Jaccard(kf, kp, nf))
      print('{} {}b {} 计算完成'.format(pid, vid, item))
  # 求max-average
  for tmp in susdic_Oc.keys():
    susdic_Oc[tmp] = susdic_Oc[tmp] - avg_susdic_Oc[tmp] / num_susdic_Oc[tmp]
  for tmp in susdic_Ds.keys():
    susdic_Ds[tmp] = susdic_Ds[tmp] - avg_susdic_Ds[tmp] / num_susdic_Ds[tmp]
  for tmp in susdic_Ja.keys():
    susdic_Ja[tmp] = susdic_Ja[tmp] - avg_susdic_Ja[tmp] / num_susdic_Ja[tmp]
  for tmp in susdic_Op.keys():
    susdic_Op[tmp] = susdic_Op[tmp] - avg_susdic_Op[tmp] / num_susdic_Op[tmp]
  for tmp in susdic_Ta.keys():
    susdic_Ta[tmp] = susdic_Ta[tmp] - avg_susdic_Ta[tmp] / num_susdic_Ta[tmp]
  for tmp in susdic_Gp.keys():
    susdic_Gp[tmp] = susdic_Gp[tmp] - avg_susdic_Gp[tmp] / num_susdic_Gp[tmp]

  # 将怀疑度列表排序
  sorted_dict_Oc = dict(sorted(susdic_Oc.items(), key=lambda item: item[1], reverse=True))
  sorted_dict_Ds = dict(sorted(susdic_Ds.items(), key=lambda item: item[1], reverse=True))
  sorted_dict_Ja = dict(sorted(susdic_Ja.items(), key=lambda item: item[1], reverse=True))
  sorted_dict_Op = dict(sorted(susdic_Op.items(), key=lambda item: item[1], reverse=True))
  sorted_dict_Ta = dict(sorted(susdic_Ta.items(), key=lambda item: item[1], reverse=True))
  sorted_dict_Gp = dict(sorted(susdic_Gp.items(), key=lambda item: item[1], reverse=True))

  # 将排序后的字典写入csv文件,选择对应的文件
  # myfile = '/home/changzexing/susMaxSubAvgMbertType1/' + pid + '/' + vid + 'b'
  # myfile = '/home/changzexing/susMaxSubAvgMbertType3/' + pid + '/' + vid + 'b'
  # myfile = '/home/changzexing/susMaxSubAvgMbertType4/' + pid + '/' + vid + 'b'
  myfile = '/home/changzexing/susMaxSFClu/' + pid + '/' + vid + 'b'
  # myfile = '/home/changzexing/susMaxSubAvgMajorType1/' + pid + '/' + vid + 'b'
  # myfile = '/home/changzexing/susMaxSubAvgMajorType3/' + pid + '/' + vid + 'b'
  # myfile = '/home/changzexing/susMaxSubAvgMajorType4/' + pid + '/' + vid + 'b'
  # myfile = '/home/changzexing/susMaxMbertType1/' + pid + '/' + vid + 'b'
  # myfile = '/home/changzexing/susMaxMbertType3/' + pid + '/' + vid + 'b'
  # myfile = '/home/changzexing/susMaxMbertType4/' + pid + '/' + vid + 'b'
  # myfile = '/home/changzexing/susMaxMajorType1/' + pid + '/' + vid + 'b'
  # myfile = '/home/changzexing/susMaxMajorType3/' + pid + '/' + vid + 'b'
  # myfile = '/home/changzexing/susMaxMajorType4/' + pid + '/' + vid + 'b'
  print(myfile)
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
     


if __name__ == '__main__':
    # 该脚本负责计算高阶变异体的怀疑度
    # 1. 确定309/321文件位置
    # 2. 统一183和257的类型，要求一样，例如183是type3.257也需要是type3
    # 3. 如果聚合方式为max-avg，则需要打开229-240行，注释时是max的情况
    # pid = 'Time'
    # svid = 1
    # evid = 27
    pid = 'Chart'
    svid = 1
    evid = 26
    # pid = 'Math'
    # svid = 1
    # evid = 106
    # pid = 'Math'
    # svid = 1
    # evid = 106
    name = "mutant_result_faulty_file_SFClu_shanjian_json"
    for i in range(svid, evid + 1):
      # if str(i) == '2':
      #     continue
      # list_ignore = [5, 6, 13, 14, 15, 16, 17, 19, 36, 44, 54, 59, 71, 73, 74, 99]
      # if i in list_ignore:
      #    continue
      # init('/home/changzexing/d4jbasecover/' + pid + '/' + str(i) + 'b', '/home/changzexing/d4jclean/' + pid + '/' + str(i) + 'b', '/home/changzexing/mutant_result_faulty_file_hebing_json/' + pid + '/' + str(i) + 'b', str(i), pid)
      init(f'/home/changzexing/d4jbasecover/{pid}/{i}b', f'/home/changzexing/d4jclean/{pid}/{i}b', f'/home/changzexing/{name}/{pid}/{i}b', str(i), pid)
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

