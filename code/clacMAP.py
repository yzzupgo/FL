import json
import os
import pandas as pd
class DataHandler:
    def __init__(self, sheet_names):
        # Initialize a DataFrame for each sheet
        self.dfs = {name: pd.DataFrame(columns=['Project',	'Granularity',	'Formula',	'Aggregation',	'MutationTool',	'TieBreak',	'Type','Map','FileName']) for name in sheet_names}

    def add_data(self, sheet_name, data_dict):
        # Append new data to the specified DataFrame
        self.dfs[sheet_name] = self.dfs[sheet_name].append(data_dict, ignore_index=True)

    def save_data(self, filename):
        # Save all DataFrames to an Excel file, each DataFrame to a different sheet
        with pd.ExcelWriter(filename) as writer:
            for name, df in self.dfs.items():
                df.to_excel(writer, sheet_name=name, index=False)


def read_json_file(file_path):
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
            return data
    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f"Error: {e}")
        return {}

def init(handler):
  # Delta4MsMajorSFCluType3RankAve MaxMajorLast2FirstStatementType3RankAve RandomMixStatementType3RankAve
  # filenameList = ["major", "mBert", "merge"]
#   approachList = ["Chart", "Time", "Lang", "Math", "Mockito", "Closure", "Cli", "Codec", "Gson", "JacksonXml", "Csv"]
  approachList = ["Chart", "Time", "Lang", "Math", "Mockito", "Closure"]
  # filenameList = ["Major", "Mbert"]
  filenameList = ["Major"]
  #formula = ["Dstar", "Ochiai", "Jaccard", "Op2", "Tarantula", "Gp13"]
  formula = ["Dstar", "Ochiai", "Gp13"]
  #aggregation_list = ["MaxMajorLast2FirstStatement", "RandomMixStatement", "Delta4MsMajorSFClu"]
 
  #aggregation_list = ["MaxMajor", "MuseMaxMajor", "MaxMajorSamelineFOM","Delta4MsMajorSFClu"]

  aggregation_list = ["MaxMajor", "MuseMaxMajor", "Delta4MsMajorSFClu"]
  # aggregation_list = ["MaxMajor", "MaxMbert", "MuseMaxMajor", "MuseMaxMbert"]
  # aggregation_list = ["Delta4MsMajorSFClu", "Delta4MsMajorMbertSFClu"]
  #aggregation_list = ["Delta4MsMajorSFClu","MaxMajorSFDis","MaxMajorSamelineFOM"]
  type_list = ["Type3"] #["Type1","Type3","Type4"]
  Tie_breakList = ["Ave"] # ["Ave","Best"]
  # formula = ["Dstar_function", "Ochiai_function", "Jaccard_function", "Op2_function", "Tarantula_function", "Gp13_function"]
  ans = {}
  
  for filename in filenameList:  # 两个方法
    ans[filename] = {}
    for item in formula: # 每个公式
      for aggregation in aggregation_list: # 每个聚合方法
        for type in type_list: # 每个type
          for Tie_break in Tie_breakList: # 每个tie_break
            mbert_map = major_map = 0
            num_mbert = num_major = 0
            for approach in approachList: # 每个项目
              print(f" {item} {aggregation} {type} {Tie_break}")
              # MaxSubAvgMajorType1RankAve 
              faulty_line_rank_path = f"/home/changzexing/{aggregation}{type}Rank{Tie_break}/{approach}/{item}.json" # 这个文件那里来的
              faulty_line_rank = read_json_file(faulty_line_rank_path) # 错误行的rank
              sum = 0
              num = 0
              for vid in faulty_line_rank: # 对于每个版本
                tmp = faulty_line_rank[vid] # 这个版本的所有的   错误行：rank
                num += len(tmp)
                for line in tmp:
                  sum += 1 / tmp[line]
                  if filename == "Mbert":
                    mbert_map += 1 / tmp[line]
                    num_mbert += 1
                  else:
                    major_map += 1 / tmp[line]
                    num_major += 1
              if num == 0:
                ans[filename][item] = 0
              else:
                ans[filename][item] = sum / num
              dict_tmp = {
                "Project": approach,
                "Granularity": "Line",
                "Formula": item,
                "Aggregation": aggregation,
                "MutationTool": filename,
                "TieBreak": Tie_break,
                "Type": type,
                "Map": ans[filename][item],
                "FileName": f"{aggregation}{type}Rank{Tie_break}"
              }
              handler.add_data("sheet1", dict_tmp)
            if filename == "Mbert":
              dict_tmp = {
                "Project": "All",
                "Granularity": "Line",
                "Formula": item,
                "Aggregation": aggregation,
                "MutationTool": filename,
                "TieBreak": Tie_break,
                "Type": type,
                "Map": mbert_map / num_mbert,
                "FileName": f"{aggregation}{type}Rank{Tie_break}"
              }
              handler.add_data("sheet1", dict_tmp)
            else:
              dict_tmp = {
                "Project": "All",
                "Granularity": "Line",
                "Formula": item,
                "Aggregation": aggregation,
                "MutationTool": filename,
                "TieBreak": Tie_break,
                "Type": type,
                "Map": major_map / num_major,
                "FileName": f"{aggregation}{type}Rank{Tie_break}"
              }
              handler.add_data("sheet1", dict_tmp)

  
  # map_json_path = f"/home/changzexing/Map{pro}"
  # if not os.path.exists(map_json_path):
  #   os.makedirs(map_json_path)
  # # json_str = json.dumps(ans)
  # with open(f"{map_json_path}/{pro}.json", 'w') as f:
  #   json.dump(ans, f)
      


if __name__ == "__main__":
    sheet_names = ["sheet1"]  # Specify your sheet names
    handler = DataHandler(sheet_names)
    init(handler)
    # handler.save_data("/home/changzexing/d4jscript/hblscript/MAP/hblclacMAPMetallaxis.xlsx")  # Specify your filename
    handler.save_data("/home/changzexing/d4jscript/hblscript/MAP/MajorDeltaMs.xlsx")