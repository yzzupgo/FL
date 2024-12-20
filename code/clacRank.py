import csv
from fileinput import filename
import os
import json

#读取txt
def read_txt_file(file_path):
    lines = []
    with open(file_path, 'r') as f:
        for line in f:
            lines.append(line.strip())
    return lines

#将错误行的信息处理成和sus列表中的格式相同
def changeTxtEqualCsv(mutant_dic):
    faultyList = []
    for key, value in mutant_dic.items():
    #   faultyLine = key[len(sourceFilePath):].split(".")[0].replace("/", "-")
      faultyLine = key[1:].split(".")[0].replace("/", "-")
      for lineNum in value:
        # 之后要重写 这是转储命名有问题
        # tmp = faultyLine + "-" + str(lineNum)
        # tmp = tmp.split('-')
        # del tmp[-2]
        # tmp = '-'.join(tmp)
        # faultyList.append(tmp)
        # 到这里结束
        faultyList.append(faultyLine + "-" + str(lineNum))
    return faultyList

#获取rank及csv数据,下面有问题
def get_rank_statement(csv_path):
  # 存储value值与排名的映射关系
  value_rank_dict = {}
  with open(csv_path, 'r') as f:
      reader = csv.reader(f)
      # next(reader) # 跳过header
      sus = {}
      rank = 1
      num = 1
      for row in reader:
          sus[row[0]] = row[1]
          value = str(row[1])
          if value not in value_rank_dict:
              # rank = num
              value_rank_dict[value] = rank
              rank += 1
          # num += 1
  return sus, value_rank_dict

#读取权重文件夹名称
def get_subdirectories(directory):
    subdirs = []
    for root, dirs, files in os.walk(directory):
        # 这里只关心根目录下的直接子目录
        if root == directory:
            subdirs.extend(dirs)
            break  # 如果你只需要遍历最顶层目录的子文件夹，可以使用break提前退出循环
    return subdirs

#生成csv
def load_csv(path, content, formula, rank, faultylist):
  try:
    ans = {} # type: dict[str, dict[str, int]]
    for key, value in content.items(): # key: "src-main-java-org-apache-commons-lang3-math-NumberUtils-1414", value: {"rank": 1, "faulty": false}
      ans[key] = {
        "rank": rank[value],
        "faulty": True if key in faultylist else False,
        "sus": value
      }
    with open(f"{path}/{formula}.json", 'w') as f:
      json.dump(ans, f)
    return True
  except:
    return False

def convert_str_to_dict(s):
    start_index = 0
    result = {}
    while True:
        key_start_index = s.find("'", start_index) + 1
        if key_start_index == 0:
            break
        key_end_index = s.find("'", key_start_index)
        key = s[key_start_index:key_end_index]
        value_start_index = s.find("[", key_end_index) + 1
        value_end_index = s.find("]", value_start_index)
        value_str = s[value_start_index:value_end_index]
        if value_str:
            value = list(map(int, value_str.split(",")))
            result[key] = value
        start_index = value_end_index
    return result
  
def init(pro, svid, evid): 
  directory_path = '/home/changzexing/weighted2_results'
  subdirectories = get_subdirectories(directory_path)
  formula = ["Dstar", "Ochiai", "Jaccard", "Op2", "Tarantula", "Gp13"]
  # formula = ["Dstar_function", "Ochiai_function", "Jaccard_function", "Op2_function", "Tarantula_function", "Gp13_function"]
  # filename = ["susNew", "susMajor", "susMerge"]
  # filename = ["susNew", "susMajor"]
  # filename = ["susMaxSubAvgMbert"]
  # filename = ["susMaxMbertType1", "susMaxMbertType3", "susMaxMbertType4"]
  # filename = ["susMaxSubAvgMbert", "susMaxSubAvgMbertType1", "susMaxSubAvgMbertType3"]
  # filename = ["susMaxSubAvgMajorType1", "susMaxSubAvgMajorType3", "susMaxSubAvgMajorType4"]
  # filename = ["susMaxMajorType1"]
  # filename = ['susMaxMajorType1','susMaxMbertType1']
  # filename = ['susMaxMajorType3','susMaxMbertType3', 'MusesusMaxMajorType3', 'MusesusMaxMbertType3']
  # filename = ['susMaxMajorLast2FirstStatementType3', 'susMaxMajorRandomMixStatementType3', 'susMaxMajorSFCluDelta4MsType3']
  # filename = ['MusesusMaxMbertType3','MusesusMaxMajorType3']
  # filename = ['susMaxMajorSFCluDelta4MsType3']
  # filename = ['susMajorMbertSFCluDelta4MsType3']
  # filename = ['susMaxMajorDelta4MsType3']
  # filename = ['MusesusMaxMajorMbertType3', 'MesusMaxMajorMbertType3']
  # filename = ['susMaxMajorMbertSamelineStatementType3']
  filename = ['susMaxMajorSamelineFOMType3', 'susMaxMajorSFDisType3']
  # filename为待转换的csv文件夹，tmp为转化后的rank的文件夹名称
  # filename = ["SusMajorMergeMaxType3","SusMbertMergeMaxType3","SusMergeMaxType3"]
  # filename = ["susMaxMbertType1", "susMaxMbertType3", "susMaxMbertType4", "susMaxSubAvgMbertType4", "susMaxSubAvgMbertType1", "susMaxSubAvgMbertType3", "susMaxSubAvgMajorType1", "susMaxSubAvgMajorType3", "susMaxSubAvgMajorType4", "susMaxMajorType1", "susMaxMajorType3", "susMaxMajorType4"]
  # 这里只需要生成一个文件夹就可以了susMaxMergeType3，susMaxMbertMergeType3，susMaxMajorMergeType3
  # filename = ["susMaxMbertType1_func", "susMaxMbertType3_func", "susMaxMbertType4_func", "susMaxSubAvgMbertType4_func", "susMaxSubAvgMbertType1_func", "susMaxSubAvgMbertType3_func", "susMaxSubAvgMajorType1_func", "susMaxSubAvgMajorType3_func", "susMaxSubAvgMajorType4_func", "susMaxMajorType1_func", "susMaxMajorType3_func", "susMaxMajorType4_func"]
  # filename = ["susMaxMbertType1", "susMaxMajorType1"] 
  # filename = ["susMaxSFCluShanjianFOM"]
  # filename = ['susMaxMajorSFCluFOM']
  # filename = ["susMaxSFCluQuanjiFOM"]
  # filename = ['susMaxMajorLast2FirstFOM', 'susMaxMajorRandomMixFOM']
  # filename = ['susMaxMajorMbertLast2FirstFOM', 'susMaxMajorMbertRandomMixFOM']
  # filename = ['susMaxMajorMbertSFCluFOM', 'susMaxMajorMbertSFCluFOM']
  tmp = {
     'susMaxMajorSFDisType3': 'susMaxMajorSFDisType3Rank',
     'susMaxMajorMbertSamelineStatementType3': 'susMaxMajorMbertSamelineStatementType3Rank',
     'susMaxMajorSamelineFOMType3': 'susMaxMajorSamelineFOMType3Rank',
     'susMaxMajorLast2FirstStatementType3': 'susMaxMajorLast2FirstStatementType3Rank', 
     'susMaxMajorRandomMixStatementType3': 'susMaxMajorRandomMixStatementType3Rank', 
     'susMaxMajorSFCluDelta4MsType3': 'susMaxMajorSFCluDelta4MsType3Rank',
     'susMaxMajorDelta4MsType3': 'susMaxMajorDelta4MsType3Rank',
     "susMajorMbertSFCluDelta4MsType3" : "susMajorMbertSFCluDelta4MsType3Rank",
     "susMaxMbertHebingSingleType1" : "susMaxMbertHebingSingleType1Rank",
     "susMaxMajorMberthebing": "susMaxMajorMberthebingRank",
     "susMaxSFCluSingle": "susMaxSFCluSingleRank",
     "susMaxSFCluShanjianFOM": "susMaxSFCluShanjianFOMRank",
     "susMaxMajorSFCluFOM": "susMaxMajorSFCluFOMRank",
     "susMaxSFCluQuanjiFOM": "susMaxSFCluQuanjiFOMRank",
     "susMaxMajorLast2FirstFOM": "susMaxMajorLast2FirstFOMRank",
     "susMaxMajorRandomMixFOM": "susMaxMajorRandomMixFOMRank",
     "susMaxMajorMbertLast2FirstFOM": "susMaxMajorMbertLast2FirstFOMRank",
     "susMaxMajorMbertRandomMixFOM": "susMaxMajorMbertRandomMixFOMRank",
     "susMaxMajorMbertSFCluFOM": "susMaxMajorMbertSFCluFOMRank",
     "susMaxMajorMbertSFCluFOMType3": "susMaxMajorMbertSFCluFOMType3Rank",

     "SusMajorMergeMaxType3" : "SusMajorMergeMaxType3Rank",
     "SusMbertMergeMaxType3" : "SusMbertMergeMaxType3Rank",
     "SusMergeMaxType3" : "SusMergeMaxType3Rank",
     "susNew": "mbertrank",
     "susMajor": "majorrank",
     "susMerge": "mergerank",
     "susMaxSubAvgMbertType4": "susMaxSubAvgMbertType4Rank",
     "susMaxSubAvgMbertType1": "susMaxSubAvgMbertType1Rank",
     "susMaxSubAvgMbertType3": "susMaxSubAvgMbertType3Rank",
     "susMaxMbertType1" : "susMaxMbertType1Rank",
     "susMaxMbertType3" : "susMaxMbertType3Rank",
     "susMaxMbertType4" : "susMaxMbertType4Rank",
     "susMaxSubAvgMajorType1" : "susMaxSubAvgMajorType1Rank",
     "susMaxSubAvgMajorType3" : "susMaxSubAvgMajorType3Rank",
     "susMaxSubAvgMajorType4" : "susMaxSubAvgMajorType4Rank",
     "susMaxMajorType1" : "susMaxMajorType1Rank",
     "susMaxMajorType3" : "susMaxMajorType3Rank",
     "susMaxMajorType4" : "susMaxMajorType4Rank",

     "susMaxSubAvgMbertType4_func": "susMaxSubAvgMbertType4Rank_func",
     "susMaxSubAvgMbertType1_func": "susMaxSubAvgMbertType1Rank_func",
     "susMaxSubAvgMbertType3_func": "susMaxSubAvgMbertType3Rank_func",
     "susMaxMbertType1_func" : "susMaxMbertType1Rank_func",
     "susMaxMbertType3_func" : "susMaxMbertType3Rank_func",
     "susMaxMbertType4_func" : "susMaxMbertType4Rank_func",
     "susMaxSubAvgMajorType1_func" : "susMaxSubAvgMajorType1Rank_func",
     "susMaxSubAvgMajorType3_func" : "susMaxSubAvgMajorType3Rank_func",
     "susMaxSubAvgMajorType4_func" : "susMaxSubAvgMajorType4Rank_func",
     "susMaxMajorType1_func" : "susMaxMajorType1Rank_func",
     "susMaxMajorType3_func" : "susMaxMajorType3Rank_func",
     "susMaxMajorType4_func" : "susMaxMajorType4Rank_func",

     "MusesusMaxMbertType3" : "MusesusMaxMbertType3Rank",
     "MusesusMaxMajorType3" : "MusesusMaxMajorType3Rank",
     "MusesusMaxMajorMbertType3": "MusesusMaxMajorMbertType3Rank",
     "MesusMaxMajorMbertType3": "MesusMaxMajorMbertType3Rank"
  }
  mutant_path = "/home/changzexing/faultyLine/{}FalutLine.txt".format(pro)# 语句级错误行的路径
  # mutant_path = "/home/changzexing/faultyFunc/{}FalutLine_function.txt".format(pro)# 方法错误行的路径
  faultyInfo = read_txt_file(mutant_path) # 错误行信息!
  faultyList = {}
  for _item in faultyInfo:
    vid = _item.split(" ")[1]
    # 获取错误行信息
    # mutant = item[(item.index('{') + 1):-1]
    mutant = _item[(_item.index('{')):]
    mutant_dic = convert_str_to_dict(mutant)
    faultyList[vid] = changeTxtEqualCsv(mutant_dic)
  for item in range(svid, evid + 1):
    list_ignore = [] 
    # if pro == 'Math':
    #     list_ignore = [5, 6, 13, 14, 15, 16, 17, 19, 36, 44, 54, 59, 71, 73, 74, 99]
    # elif pro == 'Chart':
    #     list_ignore = [8, 20] #12 和20是hebing版本加的
    # elif pro == 'Lang':
    #     list_ignore = [56, 30, 22,39,20,40,29,59,18,45,43,8,14,38,26,31]
    if item in list_ignore:
        print(f"跳过{item}版本")
        continue
    for formula_item in formula:
      for _item in filename:
        mbert_sus_path = f"/home/changzexing/{_item}/{pro}/{item}b/{formula_item}.csv" # 每个版本的怀疑度csv文件的路径,我的文件忘记加b了，记得加
        print(mbert_sus_path)
        if not os.path.exists(mbert_sus_path):
          print("不存在")
          continue
        statement, mbert_rank = get_rank_statement(mbert_sus_path) # 传入csv文件，输出
        mbert_rank_path = f"/home/changzexing/{tmp[_item]}/{pro}/{item}b" # 生成rank值的目标路径
        print(mbert_rank_path)
        if not os.path.exists(mbert_rank_path):
          os.makedirs(mbert_rank_path)
        if str(item) not in faultyList:
          continue 
        result = load_csv(mbert_rank_path, statement, formula_item, mbert_rank, faultyList[str(item)])
        print(f"{pro}-{item}b-{_item}", "完成" if result else "失败")
      # break
    # break



if __name__ == '__main__':
  init('Chart', 1, 26)
  # init('Time', 1, 27)
  # init('Lang', 1, 65)
  # init('Math', 1, 106)
  # init('Mockito', 1, 38)
  # init('Closure', 1, 133)
  # init('Cli', 1, 39)
  # init('Codec', 1, 18)
  # init('Csv', 1, 16) 
  # init('Gson', 1, 18)
  # init('JacksonXml', 1, 6)