import csv
import json
import os

#csv转字典
def csv_to_dict(file_path):
    json_content = read_json_file(file_path)
    ranknum = {}
    trueranknum = {}
    faulty = []
    sum = 0
    for item in json_content:
      tmp_rank = json_content[item]["rank"]
      tmp_rank = str(tmp_rank)
      if json_content[item]["faulty"] == False:
        if tmp_rank in trueranknum:
            trueranknum[tmp_rank].append(item)
        else:
            trueranknum[tmp_rank] = [item]
      if tmp_rank in ranknum:
          ranknum[tmp_rank].append(item)
      else:
          ranknum[tmp_rank] = [item]
      if json_content[item]["faulty"] == True:
          faulty.append(item)
    return json_content, ranknum, faulty, trueranknum, sum

def read_json_file(file_path):
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
            return data
    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f"Error: {e}")
        return {}

def init(pro, svid, evid, filename):
  formula = ["Dstar", "Ochiai", "Jaccard", "Op2", "Tarantula", "Gp13"]
  ans = {}
  faultyFileLines_path = f"/home/changzexing/faultyFileLines/{pro}.json" # 语句级的
  faultyFileLines = read_json_file(faultyFileLines_path)
  for item in range(svid, evid + 1):
    list_notrun_chart = [8,10]
    if item in list_notrun_chart:
      print(f"跳过{item}版本")
    # list_notrun_lang = [56,30,22,39,20,40,29,59,18,45,43,8,14,38,26,31]
    # if item in list_notrun_lang:
    #   print(f"跳过{item}版本")
    #   continue
    # list_ignore_math = [5, 6, 13, 14, 15, 16, 17, 19, 36, 44, 54, 59, 71, 73, 74, 99]
    # if item in list_ignore_math:
    #     print(f"跳过{item}版本")
    #     continue  
    for formula_item in formula:
      mbert_sus_path = f"/home/changzexing/{filename}rank/{pro}/{item}b/{formula_item}.json"
      if not os.path.exists(mbert_sus_path):
        continue
      sus, ranknum, faulty, trueranknum, sum = csv_to_dict(mbert_sus_path)
      for _item in faulty:
        i = 0
        sus[_item]["rank"]  = str(sus[_item]["rank"])
        if int(sus[_item]["rank"]) > 1:
          for rank in range(1, int(sus[_item]["rank"])):
            i += len(trueranknum[str(rank)]) if str(rank) in trueranknum else 0
        j = len(ranknum[sus[_item]["rank"]]) if sus[_item]["rank"] in ranknum else 0
        rank = ((i + 1) + (i + j)) / 2
        examval = rank / faultyFileLines[f"{pro}-{item}"]
        if formula_item in ans:
          ans[formula_item].append(examval)
        else:
          ans[formula_item] = [examval]
  exam_json_path = f"/home/changzexing/Exam{filename}"
  if not os.path.exists(exam_json_path):
    os.makedirs(exam_json_path)
  with open(f"{exam_json_path}/{pro}.json", 'w') as f:
    json.dump(ans, f)


if __name__ == '__main__':
  # filename = ["major", "mbert", "merge"]
  filename = ["major", "mbert"]
  for item in filename:
    # init('Math', 1, 106, item)
    init('Chart', 1, 26, item)
    # init('Time', 1, 27, item)
    # init('Lang', 1, 65, item)