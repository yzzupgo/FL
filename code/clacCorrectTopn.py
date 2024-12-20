import json
import os

def read_json_file(file_path):
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
            return data
    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f"Error: {e}")
        return {}

def init(pro):


    # 定义每个项目类型的有效版本号
  # valid_versions = {
  #       "Chart": ["18b", "17b", "5b", "7b"],
  #       "Time": ["13b", "1b", "23b"],
  #       "Closure": ["100b", "105b", "110b", "21b", "25b", "32b", "34b", "39b", "45b", "49b", "4b", "50b", "54b", "68b", "6b", "72b", "75b", "76b", "85b", "90b", "99b", "9b", "14b", "115b", "131b"],
  #       "Lang": ["15b", "19b", "20b", "27b", "32b", "34b", "35b", "36b", "41b", "4b", "50b", "56b", "60b", "7b", "40b", "42b", "46b", "25b"],
  #       "Math": ["100b", "14b", "18b", "24b", "26b", "29b", "35b", "37b", "38b", "44b", "46b", "49b", "54b", "62b", "64b", "65b", "66b", "68b", "6b", "72b", "76b", "7b", "81b", "83b", "86b", "92b", "93b", "55b", "56b", "57b", "58b", "59b", "5b", "60b", "43b"],
  #       "Mockito": ["4b", "11b"]
  #   }
#   valid_versions = {
#     "Chart": ['16b', '22b', '25b', '1b', '6b', '8b'],
#     "Time": ['12b', '2b', '10b', '11b', '17b', '19b', '20b', '22b', '5b', '7b', '8b', '26b', '9b'],
#     "Closure": ['106b', '30b', '37b', '47b', '79b', '89b', '16b', '17b', '18b', '20b', '22b', '23b', '24b', '27b', '31b', '35b', '38b', '3b', '40b', '46b', '48b', '51b', '52b', '55b', '57b', '59b', '62b', '64b', '65b', '67b', '70b', '71b', '108b', '112b', '117b', '121b', '122b', '126b', '127b', '128b', '7b', '87b', '97b'],
#     "Lang": ['30b', '16b', '17b', '53b', '57b', '58b', '59b', '61b', '63b', '6b', '8b', '10b', '18b', '1b', '22b'],
#     "Math": ['22b', '98b', '61b', '63b', '67b', '69b', '70b', '74b', '75b', '77b', '79b', '80b', '82b', '85b', '87b', '88b', '8b', '91b', '94b', '9b', '97b', '15b', '16b', '21b', '23b', '31b', '40b', '42b', '47b', '52b'],
#     "Mockito": ['10b', '14b', '16b', '17b', '23b', '25b', '32b', '35b', '6b', '13b', '19b', '26b', '27b', '28b', '29b', '30b', '33b', '34b', '38b', '5b', '8b', '1b', '20b', '21b', '3b']
# }
  valid_versions = {
    "Chart": ['13b', '20b', '24b'],
    "Time": ['16b', '4b'],
    "Closure": ['73b', '74b', '78b', '109b', '111b', '114b', '125b', '130b', '132b', '83b', '82b', '96b', '101b'],
    "Lang": [],
    "Math": ['95b', '96b', '101b', '105b', '12b', '20b', '2b', '41b', '50b'],
    "Mockito": ['22b', '24b']
}



    # 根据项目类型获取有效版本
  versions_to_load = valid_versions.get(pro, [])
  # filename = ["major", "mBert", "merge"]
  # filename = ["mBert"]
  formula = ["Dstar", "Ochiai", "Jaccard", "Op2", "Tarantula", "Gp13"]
  # formula = ["Dstar_function", "Ochiai_function", "Jaccard_function", "Op2_function", "Tarantula_function", "Gp13_function"]
  # averank = "MaxSubAvgAveRang"
  # avetopn = "MaxSubAvgAveTopn"
  # bestrank = "MaxSubAvgBestRang"
  # besttopn = "MaxSubAvgBestTopn"
  
  #语句级的怀疑度：
  averank = {
    # "0": "MaxSubAvgMbertType4RankAve",
    # "1": "MaxSubAvgMbertType1RankAve",
    # "2": "MaxSubAvgMbertType3RankAve",
    # "3" : "MaxMbertType1RankAve",
    # "4" : "MaxMbertType3RankAve",
    # "5" : "MaxMbertType4RankAve",
    # "6" : "MaxSubAvgMajorType1RankAve",
    # "7" : "MaxSubAvgMajorType3RankAve",
    # "8" : "MaxSubAvgMajorType4RankAve",
    # "9" : "MaxMajorType1RankAve",
    # "10" : "MaxMajorType3RankAve",
    # "11" : "MaxMajorType4RankAve",
    # "12" : "MaxMajorMergeType3RankAve",
    # "13":  "MaxMbertMergeType3RankAve",
    # "14":  "MaxMergeType3RankAve",
    # "15": "MaxMbertHebingSingleType1RankAve",
    # "16": "MaxMajorMberthebingRankAve",
    # "17": "MaxSFCluSingleRankAve",
    # "18": 'MaxSFCluShanjianFOMRankAve',
    # "19":  "MaxMajorLast2FirstFOMRankAve",
    # "20": "MaxMajorRandomMixFOMRankAve",
    # "21": "MaxMajorSFCluFOMRankAveNew",
    # "22": "MaxMajorMbertLast2FirstFOMRankAve",
    # "23": "MaxMajorMbertRandomMixFOMRankAve",
    # "24": "MaxMajorMbertSFCluFOMRankAveNew",
    #  "25": "MaxMajorMbertSFCluFOMType3RankAveNew",
    #  "26": "Delta4MsMajorMbertSFCluType3RankAve",
    #  "27" : "MuseMaxMbertType3RankAve",
    # "28" : "MuseMaxMajorType3RankAve",
    #  "29": "MaxMajorLast2FirstStatementType3RankAve",
    #  "30": "RandomMixStatementType3RankAve",
      "31": "Delta4MsMajorSFCluType3RankAve",
     "32": "MaxMajorSamelineFOMType3RankAve",
    # "33": "MuseMaxMajorMbertType3RankAve",
    # "34": 'MaxMajorMbertSamelineStatementType3RankAve',
    # "35": 'MeMaxMajorMbertType3RankAve',
    # "36": "MaxMajorDelta4MsType3RankAve",
    "37": "MaxMajorSFDisType3RankAve"

  }
  avetopn = {
    "0": "MaxSubAvgMbertType4TopnAve",
    "1": "MaxSubAvgMbertType1TopnAve",
    "2": "MaxSubAvgMbertType3TopnAve",
    "3" : "MaxMbertType1TopnAve",
    "4" : "MaxMbertType3TopnAve",
    "5" : "MaxMbertType4TopnAve",
    "6" : "MaxSubAvgMajorType1TopnAve",
    "7" : "MaxSubAvgMajorType3TopnAve",
    "8" : "MaxSubAvgMajorType4TopnAve",
    "9" : "MaxMajorType1TopnAve",
    "10" : "MaxMajorType3TopnAve",
    "11" : "MaxMajorType4TopnAve",
    "12" : "MaxMajorMergeType3TopnAve",
    "13":  "MaxMbertMergeType3TopnAve",
    "14":  "MaxMergeType3TopnAve",
    "15": "MaxMbertHebingSingleType1TopnAve",
    "16": "MaxMajorMberthebingTopnAve",
    "17": "MaxSFCluSingleTopnAve",
    "18": 'MaxSFCluShanjianFOMTopnAve',
    "19":  "MaxMajorLast2FirstFOMTopnAve",
    "20": "MaxMajorRandomMixFOMTopnAve",
    "21": "MaxMajorSFCluFOMTopnAveNew",
    "22": "MaxMajorMbertLast2FirstFOMTopnAve",
    "23": "MaxMajorMbertRandomMixFOMTopnAve",
    "24": "MaxMajorMbertSFCluFOMTopnAveNew",
    "25": "MaxMajorMbertSFCluFOMType3TopnAveNew",
    "26": "Delta4MsMajorMbertSFCluType3TopnAve",
    "27" : "MuseMaxMbertType3TopnAve",
    "28" : "MuseMaxMajorType3TopnAve",
    "29": "MaxMajorLast2FirstStatementType3TopnAve",
    "30": "RandomMixStatementType3TopnAve",
    "31": "Delta4MsMajorSFCluType3TopnAve",
    "32": "MaxMajorSamelineFOMType3TopnAve",
    "33": "MuseMaxMajorMbertType3TopnAve",
    "34": 'MaxMajorMbertSamelineStatementType3TopnAve',
    "35": 'MeMaxMajorMbertType3TopnAve',
    "36": "MaxMajorDelta4MsType3TopnAve",
    "37": "MaxMajorSFDisType3TopnAve"

  }
  bestrank = {
    # "0": "MaxSubAvgMbertType4RankBest",
    # "1": "MaxSubAvgMbertType1RankBest",
    # "2": "MaxSubAvgMbertType3RankBest",
    # "3" : "MaxMbertType1RankBest",
    # "4" : "MaxMbertType3RankBest",
    # "5" : "MaxMbertType4RankBest",
    # "6" : "MaxSubAvgMajorType1RankBest",
    # "7" : "MaxSubAvgMajorType3RankBest",
    # "8" : "MaxSubAvgMajorType4RankBest",
    # "9" : "MaxMajorType1RankBest",
    # "10" : "MaxMajorType3RankBest",
    # "11" : "MaxMajorType4RankBest",
    # "12" : "MaxMajorMergeType3RankBest",
    # "13":  "MaxMbertMergeType3RankBest",
    # "14":  "MaxMergeType3RankBest",
    # "15": "MaxMbertHebingSingleType1RankBest",
    # "16": "MaxMajorMberthebingRankBest",
    # "17": "MaxSFCluSingleRankBest",
    # "18": 'MaxSFCluShanjianFOMRankBest',
    "19":  "MaxMajorLast2FirstFOMRankAve",
    "20": "MaxMajorRandomMixFOMRankAve"
  }
  besttopn = {
    "0": "MaxSubAvgMbertType4TopnBest",
    "1": "MaxSubAvgMbertType1TopnBest",
    "2": "MaxSubAvgMbertType3TopnBest",
    "3" : "MaxMbertType1TopnBest",
    "4" : "MaxMbertType3TopnBest",
    "5" : "MaxMbertType4TopnBest",
    "6" : "MaxSubAvgMajorType1TopnBest",
    "7" : "MaxSubAvgMajorType3TopnBest",
    "8" : "MaxSubAvgMajorType4TopnBest",
    "9" : "MaxMajorType1TopnBest",
    "10" : "MaxMajorType3TopnBest",
    "11" : "MaxMajorType4TopnBest",
    "12" : "MaxMajorMergeType3TopnBest",
    "13":  "MaxMbertMergeType3TopnBest",
    "14":  "MaxMergeType3TopnBest",
    "15": "MaxMbertHebingSingleType1TopnBest",
    "16": "MaxMajorMberthebingTopnBest",
    "17": "MaxSFCluSingleTopnBest",
    "18": 'MaxSFCluShanjianFOMTopnBest'
  }

  #方法级的怀疑度
  # averank = {
  #   "0": "MaxSubAvgMbertType4RankAve_func",
  #   "1": "MaxSubAvgMbertType1RankAve_func",
  #   "2": "MaxSubAvgMbertType3RankAve_func",
  #   "3" : "MaxMbertType1RankAve_func",
  #   "4" : "MaxMbertType3RankAve_func",
  #   "5" : "MaxMbertType4RankAve_func",
  #   "6" : "MaxSubAvgMajorType1RankAve_func",
  #   "7" : "MaxSubAvgMajorType3RankAve_func",
  #   "8" : "MaxSubAvgMajorType4RankAve_func",
  #   "9" : "MaxMajorType1RankAve_func",
  #   "10" : "MaxMajorType3RankAve_func",
  #   "11" : "MaxMajorType4RankAve_func"
  # }
  # avetopn = {
  #   "0": "MaxSubAvgMbertType4TopnAve_func",
  #   "1": "MaxSubAvgMbertType1TopnAve_func",
  #   "2": "MaxSubAvgMbertType3TopnAve_func",
  #   "3" : "MaxMbertType1TopnAve_func",
  #   "4" : "MaxMbertType3TopnAve_func",
  #   "5" : "MaxMbertType4TopnAve_func",
  #   "6" : "MaxSubAvgMajorType1TopnAve_func",
  #   "7" : "MaxSubAvgMajorType3TopnAve_func",
  #   "8" : "MaxSubAvgMajorType4TopnAve_func",
  #   "9" : "MaxMajorType1TopnAve_func",
  #   "10" : "MaxMajorType3TopnAve_func",
  #   "11" : "MaxMajorType4TopnAve_func"
  # }
  # bestrank = {
  #   "0": "MaxSubAvgMbertType4RankBest_func",
  #   "1": "MaxSubAvgMbertType1RankBest_func",
  #   "2": "MaxSubAvgMbertType3RankBest_func",
  #   "3" : "MaxMbertType1RankBest_func",
  #   "4" : "MaxMbertType3RankBest_func",
  #   "5" : "MaxMbertType4RankBest_func",
  #   "6" : "MaxSubAvgMajorType1RankBest_func",
  #   "7" : "MaxSubAvgMajorType3RankBest_func",
  #   "8" : "MaxSubAvgMajorType4RankBest_func",
  #   "9" : "MaxMajorType1RankBest_func",
  #   "10" : "MaxMajorType3RankBest_func",
  #   "11" : "MaxMajorType4RankBest_func"
  # }
  # besttopn = {
  #   "0": "MaxSubAvgMbertType4TopnBest_func",
  #   "1": "MaxSubAvgMbertType1TopnBest_func",
  #   "2": "MaxSubAvgMbertType3TopnBest_func",
  #   "3" : "MaxMbertType1TopnBest_func",
  #   "4" : "MaxMbertType3TopnBest_func",
  #   "5" : "MaxMbertType4TopnBest_func",
  #   "6" : "MaxSubAvgMajorType1TopnBest_func",
  #   "7" : "MaxSubAvgMajorType3TopnBest_func",
  #   "8" : "MaxSubAvgMajorType4TopnBest_func",
  #   "9" : "MaxMajorType1TopnBest_func",
  #   "10" : "MaxMajorType3TopnBest_func",
  #   "11" : "MaxMajorType4TopnBest_func"
  # }
  # 计算averank
  for item in averank:
    for formula_item in formula:
      ans = {
        "Top-1": 0,
        "Top-3": 0,
        "Top-5": 0,
        "Top-10": 0,
      }
      ave_json_path = f"/home/changzexing/{averank[item]}/{pro}/{formula_item}.json"
      result = read_json_file(ave_json_path)
    #   if result:
    #             # 筛选出在有效版本列表中的键
    #             result = {k: result[k] for k in result if k in versions_to_load}
      for vid in result:
        vid = result[vid]
        for faultyline in vid:
          if vid[faultyline] <= 1:
            ans["Top-1"] += 1
          if vid[faultyline] <= 3:
            ans["Top-3"] += 1
          if vid[faultyline] <= 5:
            ans["Top-5"] += 1
          if vid[faultyline] <= 10:
            ans["Top-10"] += 1
      print(ans)
      topn_json_path = f"/home/changzexing/{avetopn[item]}/{pro}"
      if not os.path.exists(topn_json_path):
        os.makedirs(topn_json_path)
      with open(f"{topn_json_path}/{formula_item}.json", 'w') as f:
        json.dump(ans, f)
  # 计算bestrank
  # for item in bestrank:
  #   for formula_item in formula:
  #     ans = {
  #       "Top-1": 0,
  #       "Top-3": 0,
  #       "Top-5": 0,
  #       "Top-10": 0,
  #     }
  #     ave_json_path = f"/home/changzexing/{bestrank[item]}/{pro}/{formula_item}.json"
  #     result = read_json_file(ave_json_path)
  #     for vid in result:
  #       vid = result[vid]
  #       for faultyline in vid:
  #         if vid[faultyline] <= 1:
  #           ans["Top-1"] += 1
  #         if vid[faultyline] <= 3:
  #           ans["Top-3"] += 1
  #         if vid[faultyline] <= 5:
  #           ans["Top-5"] += 1
  #         if vid[faultyline] <= 10:
  #           ans["Top-10"] += 1
  #     topn_json_path = f"/home/changzexing/{besttopn[item]}/{pro}"
  #     if not os.path.exists(topn_json_path):
  #       os.makedirs(topn_json_path)
  #     with open(f"{topn_json_path}/{formula_item}.json", 'w') as f:
  #       json.dump(ans, f)
            

if __name__ == "__main__":
  init("Chart")
#   init("Time")
#   init("Lang")
#   init('Math')
#   init('Mockito')
#   init('Closure')
  # init('Cli')
  # init('Codec')
  # init('Csv') 
  # init('Gson')
  # init('JacksonXml')
