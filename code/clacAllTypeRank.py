import csv
from fileinput import filename
import json
import os
from re import T

def read_json_file(file_path):
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
            return data
    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f"Error: {e}")
        return {}
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
    # print(ranknum)
    return json_content, ranknum, faulty, trueranknum, sum

def init(pro, svid, evid):
  formula = ["Dstar", "Ochiai", "Jaccard", "Op2", "Tarantula", "Gp13"]
  # formula = ["Dstar_function", "Ochiai_function", "Jaccard_function", "Op2_function", "Tarantula_function", "Gp13_function"]
  # filename = ["susMaxSFCluQuanjiFOMRank"]
  # filename = ['susMaxMajorLast2FirstFOMRank', 'susMaxMajorRandomMixFOMRank']
  # filename = ['susMaxMajorMbertLast2FirstFOMRank', 'susMaxMajorMbertRandomMixFOMRank']
  # filename = ['susMaxMajorMbertSFCluFOMRank']
  # filename = ['susMaxMbertType1Rank', 'susMaxMajorType1Rank']
  # filename = ['susMaxMbertType3Rank', 'susMaxMajorType3Rank', 'MusesusMaxMbertType3Rank','MusesusMaxMajorType3Rank']
  # filename = ['susMaxMajorSFCluDelta4MsType3Rank', 'susMaxMajorRandomMixStatementType3Rank', 'susMaxMajorLast2FirstStatementType3Rank']
  # filename = ['MusesusMaxMbertType3Rank','MusesusMaxMajorType3Rank']
  # filename = ['susMaxMajorSFCluDelta4MsType3Rank']
  # filename = ['susMajorMbertSFCluDelta4MsType3Rank']
  # filename = ['susMaxMajorDelta4MsType3Rank']
  # filename = ['MusesusMaxMajorMbertType3Rank', 'MesusMaxMajorMbertType3Rank']
  # filename = ['susMaxMajorMbertSamelineStatementType3Rank']
  filename = ['susMaxMajorSamelineFOMType3Rank', 'susMaxMajorSFDisType3Rank']
  # filename = ["susMaxMajorSFCluFOMRank"]
  # filename = ["mbertrank", "majorrank", "mergerank"]
  # filename = ["MaxSubAvgMbertrank"]
  # filename = ["susMaxSubAvgMbertType1Rank", "susMaxSubAvgMbertType3Rank"]
  # filename = ["SusMajorMergeMaxType3Rank","SusMbertMergeMaxType3Rank","SusMergeMaxType3Rank"]
  # filename = ["susMaxMbertType1Rank", "susMaxMbertType3Rank", "susMaxMbertType4Rank", "susMaxSubAvgMbertType4Rank", "susMaxSubAvgMbertType1Rank", "susMaxSubAvgMbertType3Rank", "susMaxSubAvgMajorType1Rank", "susMaxSubAvgMajorType3Rank", "susMaxSubAvgMajorType4Rank", "susMaxMajorType1Rank", "susMaxMajorType3Rank", "susMaxMajorType4Rank"]
  # filename = ["susMaxMbertType1Rank_func", "susMaxMbertType3Rank_func", "susMaxMbertType4Rank_func", "susMaxSubAvgMbertType4Rank_func", "susMaxSubAvgMbertType1Rank_func", "susMaxSubAvgMbertType3Rank_func", "susMaxSubAvgMajorType1Rank_func", "susMaxSubAvgMajorType3Rank_func", "susMaxSubAvgMajorType4Rank_func", "susMaxMajorType1Rank_func", "susMaxMajorType3Rank_func", "susMaxMajorType4Rank_func"]
  # tmp = {
  #   "mbertrank": "mBert", 
  #   "majorrank": "major", 
  #   "mergerank": "merge",
  #   "susMaxSubAvgMbertType4Rank": "mbert",
  #   "susMaxSubAvgMbertType1Rank": "mbert",
  #   "susMaxSubAvgMbertType3Rank": "mbert",
  #   "susMaxMbertType1Rank" : "mbert",
  #   "susMaxMbertType3Rank" : "mbert",
  #   "susMaxMbertType4Rank" : "mbert",
  #   "susMaxSubAvgMajorType1Rank" : "major",
  #   "susMaxSubAvgMajorType3Rank" : "major",
  #   "susMaxSubAvgMajorType4Rank" : "major",
  #   "susMaxMajorType1Rank" : "major",
  #   "susMaxMajorType3Rank" : "major",
  #   "susMaxMajorType4Rank" : "major"
  # }
  ave = {
     'susMaxMajorDelta4MsType3Rank': 'MaxMajorDelta4MsType3RankAve',
     'susMaxMajorMbertSamelineStatementType3Rank': 'MaxMajorMbertSamelineStatementType3RankAve',
     "susMaxMajorSamelineFOMType3Rank": "MaxMajorSamelineFOMType3RankAve",
     "susMaxMajorSFDisType3Rank": "MaxMajorSFDisType3RankAve",
     "susMaxMajorLast2FirstStatementType3Rank": "MaxMajorLast2FirstStatementType3RankAve",
     "susMaxMajorRandomMixStatementType3Rank": "RandomMixStatementType3RankAve",
     "susMaxMajorSFCluDelta4MsType3Rank": "Delta4MsMajorSFCluType3RankAve",
     "susMajorMbertSFCluDelta4MsType3Rank": "Delta4MsMajorMbertSFCluType3RankAve",
     "susMaxMbertHebingSingleType1Rank" : "MaxMbertHebingSingleType1RankAve",
     "susMaxMajorMberthebingRank": "MaxMajorMberthebingRankAve",
     "susMaxSFCluSingleRank": "MaxSFCluSingleRankAve",
     "susMaxSFCluShanjianFOMRank": "MaxSFCluShanjianFOMRankAve",
     "susMaxMajorSFCluFOMRank": "MaxMajorSFCluFOMRankAve",
     "susMaxSFCluQuanjiFOMRank": "MaxSFCluQuanjiFOMRankAve",
     'susMaxMajorLast2FirstFOMRank': "MaxMajorLast2FirstFOMRankAve",
     "susMaxMajorRandomMixFOMRank": "MaxMajorRandomMixFOMRankAve",
     'susMaxMajorMbertLast2FirstFOMRank': 'MaxMajorMbertLast2FirstFOMRankAve',
     'susMaxMajorMbertRandomMixFOMRank': 'MaxMajorMbertRandomMixFOMRankAve',
     "susMaxMajorMbertSFCluFOMRank": 'MaxMajorMbertSFCluFOMRankAve',
     "susMaxMajorMbertSFCluFOMType3Rank": "MaxMajorMbertSFCluFOMType3RankAve",

    "SusMajorMergeMaxType3Rank" : "MaxMajorMergeType3RankAve",
    "SusMbertMergeMaxType3Rank" : "MaxMbertMergeType3RankAve",
    "SusMergeMaxType3Rank" :      "MaxMergeType3RankAve",
    "susMaxSubAvgMbertType4Rank": "MaxSubAvgMbertType4RankAve",
    "susMaxSubAvgMbertType1Rank": "MaxSubAvgMbertType1RankAve",
    "susMaxSubAvgMbertType3Rank": "MaxSubAvgMbertType3RankAve",
    "susMaxMbertType1Rank" : "MaxMbertType1RankAve",
    "susMaxMbertType3Rank" : "MaxMbertType3RankAve",
    "susMaxMbertType4Rank" : "MaxMbertType4RankAve",
    "susMaxSubAvgMajorType1Rank" : "MaxSubAvgMajorType1RankAve",
    "susMaxSubAvgMajorType3Rank" : "MaxSubAvgMajorType3RankAve",
    "susMaxSubAvgMajorType4Rank" : "MaxSubAvgMajorType4RankAve",
    "susMaxMajorType1Rank" : "MaxMajorType1RankAve",
    "susMaxMajorType3Rank" : "MaxMajorType3RankAve",
    "susMaxMajorType4Rank" : "MaxMajorType4RankAve",
    # Muse下面的
    "MusesusMaxMbertType3Rank" : "MuseMaxMbertType3RankAve",
    "MusesusMaxMajorType3Rank" : "MuseMaxMajorType3RankAve",
    "MusesusMaxMajorMbertType3Rank": "MuseMaxMajorMbertType3RankAve",
    "MesusMaxMajorMbertType3Rank": "MeMaxMajorMbertType3RankAve",


    "susMaxSubAvgMbertType4Rank_func": "MaxSubAvgMbertType4RankAve_func",
    "susMaxSubAvgMbertType1Rank_func": "MaxSubAvgMbertType1RankAve_func",
    "susMaxSubAvgMbertType3Rank_func": "MaxSubAvgMbertType3RankAve_func",
    "susMaxMbertType1Rank_func" : "MaxMbertType1RankAve_func",
    "susMaxMbertType3Rank_func" : "MaxMbertType3RankAve_func",
    "susMaxMbertType4Rank_func" : "MaxMbertType4RankAve_func",
    "susMaxSubAvgMajorType1Rank_func" : "MaxSubAvgMajorType1RankAve_func",
    "susMaxSubAvgMajorType3Rank_func" : "MaxSubAvgMajorType3RankAve_func",
    "susMaxSubAvgMajorType4Rank_func" : "MaxSubAvgMajorType4RankAve_func",
    "susMaxMajorType1Rank_func" : "MaxMajorType1RankAve_func",
    "susMaxMajorType3Rank_func" : "MaxMajorType3RankAve_func",
    "susMaxMajorType4Rank_func" : "MaxMajorType4RankAve_func"
  }
  best = {
     "susMajorMbertSFCluDelta4MsType3Rank": "Delta4MsMajorMbertSFCluType3RankBest",
     "susMaxMbertHebingSingleType1Rank" : "MaxMbertHebingSingleType1RankBest",
     "susMaxMajorMberthebingRank": "MaxMajorMberthebingRankBest",
     "susMaxSFCluSingleRank": "MaxSFCluSingleRankBest",
     "susMaxMajorSFCluRank": "MaxMajorSFCluRankBest",
     "susMaxSFCluShanjianFOMRank": "MaxSFCluShanjianFOMRankBest",
     "susMaxMajorSFCluFOMRank": "MaxMajorSFCluFOMRankBest",
     "susMaxSFCluQuanjiFOMRank": "MaxSFCluQuanjiFOMRankBest",
     'susMaxMajorLast2FirstFOMRank': "MaxMajorLast2FirstFOMRankBest",
     "susMaxMajorRandomMixFOMRank": "MaxMajorRandomMixFOMRankBest",
     'susMaxMajorMbertLast2FirstFOMRank': 'MaxMajorMbertLast2FirstFOMRankBest',
     'susMaxMajorMbertRandomMixFOMRank': 'MaxMajorMbertRandomMixFOMRankBest',
     "susMaxMajorMbertSFCluFOMRank": 'MaxMajorMbertSFCluFOMRankBest',
     "susMaxMajorMbertSFCluFOMType3Rank": "MaxMajorMbertSFCluFOMType3RankBest",

    "SusMajorMergeMaxType3Rank" : "MaxMajorMergeType3RankBest",
    "SusMbertMergeMaxType3Rank" : "MaxMbertMergeType3RankBest",
    "SusMergeMaxType3Rank" :      "MaxMergeType3RankBest",
    "susMaxSubAvgMbertType4Rank": "MaxSubAvgMbertType4RankBest",
    "susMaxSubAvgMbertType1Rank": "MaxSubAvgMbertType1RankBest",
    "susMaxSubAvgMbertType3Rank": "MaxSubAvgMbertType3RankBest",
    "susMaxMbertType1Rank" : "MaxMbertType1RankBest",
    "susMaxMbertType3Rank" : "MaxMbertType3RankBest",
    "susMaxMbertType4Rank" : "MaxMbertType4RankBest",
    "susMaxSubAvgMajorType1Rank" : "MaxSubAvgMajorType1RankBest",
    "susMaxSubAvgMajorType3Rank" : "MaxSubAvgMajorType3RankBest",
    "susMaxSubAvgMajorType4Rank" : "MaxSubAvgMajorType4RankBest",
    "susMaxMajorType1Rank" : "MaxMajorType1RankBest",
    "susMaxMajorType3Rank" : "MaxMajorType3RankBest",
    "susMaxMajorType4Rank" : "MaxMajorType4RankBest",
    # Muse下面的
    "MusesusMaxMbertType3Rank" : "MuseMaxMbertType3RankBest",
    "MusesusMaxMajorType3Rank" : "MuseMaxMajorType3RankBest",


    "susMaxSubAvgMbertType4Rank_func": "MaxSubAvgMbertType4RankBest_func",
    "susMaxSubAvgMbertType1Rank_func": "MaxSubAvgMbertType1RankBest_func",
    "susMaxSubAvgMbertType3Rank_func": "MaxSubAvgMbertType3RankBest_func",
    "susMaxMbertType1Rank_func" : "MaxMbertType1RankBest_func",
    "susMaxMbertType3Rank_func" : "MaxMbertType3RankBest_func",
    "susMaxMbertType4Rank_func" : "MaxMbertType4RankBest_func",
    "susMaxSubAvgMajorType1Rank_func" : "MaxSubAvgMajorType1RankBest_func",
    "susMaxSubAvgMajorType3Rank_func" : "MaxSubAvgMajorType3RankBest_func",
    "susMaxSubAvgMajorType4Rank_func" : "MaxSubAvgMajorType4RankBest_func",
    "susMaxMajorType1Rank_func" : "MaxMajorType1RankBest_func",
    "susMaxMajorType3Rank_func" : "MaxMajorType3RankBest_func",
    "susMaxMajorType4Rank_func" : "MaxMajorType4RankBest_func"
  }
  # ave = "MaxSubAvgAveRangType3"
  # best = "MaxSubAvgBestRangType3"
  ans_ave = {}
  ans_best = {}
  for _item_ in filename:
    for formula_item in formula:
      for item in range(svid, evid + 1):
        list_ignore = [] 
        if item in list_ignore:
            print(f"跳过{item}版本")
            continue
        mbert_sus_path = f"/home/changzexing/{_item_}/{pro}/{item}b/{formula_item}.json"
        if not os.path.exists(mbert_sus_path):
          print(f"跳过{pro}项目{item}版本")
          continue
        sus, ranknum, faulty, trueranknum, sum = csv_to_dict(mbert_sus_path)
        # fautly_result_best = {}
        fautly_result_ave = {}
        # faulty文件全为空
        for _item in faulty:
          i = 0
          sus[_item]["rank"]  = str(sus[_item]["rank"])
          if int(sus[_item]["rank"]) > 1:
            for rank in range(1, int(sus[_item]["rank"])):
              i += len(trueranknum[str(rank)]) if str(rank) in trueranknum else 0
          # print(sus[_item]["rank"], type(sus[_item]["rank"]))
          j = len(ranknum[sus[_item]["rank"]]) if sus[_item]["rank"] in ranknum else 0
          # 平均排名
          fautly_result_ave[_item] = ((i + 1) + (i + j)) / 2
          # best
          # fautly_result_best[_item] = i + 1
        ans_ave[f"{item}b"] = fautly_result_ave
        # ans_best[f"{item}b"] = fautly_result_best
      #   break
      # break
      mbert_averank_json_path = f"/home/changzexing/{ave[_item_]}/{pro}"
      if not os.path.exists(mbert_averank_json_path):
        os.makedirs(mbert_averank_json_path)
      with open(f"{mbert_averank_json_path}/{formula_item}.json", 'w') as f:
        json.dump(ans_ave, f)
      # mbert_bestrank_json_path = f"/home/changzexing/{best[_item_]}/{pro}"
      # if not os.path.exists(mbert_bestrank_json_path):
      #   os.makedirs(mbert_bestrank_json_path)
      # with open(f"{mbert_bestrank_json_path}/{formula_item}.json", 'w') as f:
      #   json.dump(ans_best, f)




if __name__ == '__main__':
  init('Chart', 1, 26)
  init('Time', 1, 27)
  init('Lang', 1, 65)
  init('Math', 1, 106)
  init('Mockito', 1, 38)
  init('Closure', 1, 133)
  # init('Cli', 1, 39)
  # init('Codec', 1, 18)
  # init('Csv', 1, 16) 
  # init('Gson', 1, 18)
  # init('JacksonXml', 1, 6)