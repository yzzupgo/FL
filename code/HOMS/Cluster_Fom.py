from concurrent.futures import ProcessPoolExecutor
import datetime
from importlib.resources import read_binary
import json
import os
import random
import math
import sys
# from sklearn.cluster import KMeans as sklearn_KMeans
# from sklearn.cluster import AgglomerativeClustering as sklearn_AgglomerativeClustering
# import matplotlib.pyplot as plt

# from data_codeflaws import MutationRule

# 高阶变异体生成策略
class HomGenerationStrategy:

    # 聚类中，使用同一类簇中的变异体进行类RandomMix生成
    def clusterSameClassRandomMix(self, clusterdata, fom_list):
        som_list = []

        for label in clusterdata:
            classfomlist = clusterdata[label]
            while len(classfomlist) > 0:
                fom1 = fom_list[classfomlist.pop(0)]
                for loc in random.sample(classfomlist, len(classfomlist)):
                    fom2 = fom_list[loc]
                    if not fom1['message'][0][0] == fom2['message'][0][0]:
                        if fom1['message'][0][0] < fom2['message'][0][0]:
                            som_list.append([fom1['message'][0], fom2['message'][0]])
                        else:
                            som_list.append([fom2['message'][0], fom1['message'][0]])
                        classfomlist.pop(classfomlist.index(loc))
                        break
        return som_list

    # 聚类中，使用同一类簇中的变异体进行类RandomMix生成
    def clusterDifClassRandomMix(self, clusterdata, fom_list):
        som_list = []
        for i, label_i in enumerate(clusterdata):
            for j, label_j in enumerate(clusterdata):
                if i >= j:
                    continue
                # 取上三角类簇

                classfomlist_i = clusterdata[label_i]
                classfomlist_j = clusterdata[label_j]
                for fom_i_loc in classfomlist_i:
                    fom_i = fom_list[fom_i_loc]
                    for fom_j_loc in classfomlist_j:
                        fom_j = fom_list[fom_j_loc]
                        if not MutationRule(fom_i['message'][0], fom_j['message'][0]):
                            continue
                        som = sorted([fom_i['message'][0], fom_j['message'][0]], key=lambda x: x[0])
                        som_list.append(som)

        return som_list

    # 取前50%序列fom用于生成全集
    def DistancePowerGeneration(self, clusterdata, fom_list, k=0.5):
        som_list = []
        for i, label_i in enumerate(clusterdata):
            for j, label_j in enumerate(clusterdata):
                if i > j:
                    continue
                classfomlist_i = clusterdata[label_i]
                classfomlist_j = clusterdata[label_j]
                for cut_i, fom_i_loc in enumerate(classfomlist_i):
                    if cut_i >= math.ceil(len(classfomlist_i)*k):
                        continue
                    fom_i = fom_list[fom_i_loc]
                    for cut_j, fom_j_loc in enumerate(classfomlist_j):
                        if cut_j >= math.ceil(len(classfomlist_j)*k):
                            continue
                        fom_j = fom_list[fom_j_loc]
                        if not fom_i['message'][0][0] == fom_j['message'][0][0]:
                            som_list.append(sorted([fom_i['message'][0], fom_j['message'][0]], key=lambda x: x[0]))
        return som_list





# 输入多种特征，组合输出所需特诊
class Feature:

    # 以ps为空间中心，行为差异为相应特征
    def or_center_behavior(self, or_out, mut_out):
        feature_list = []
        for i in range(len(or_out)):
            if or_out[i] == mut_out[i]:
                feature_list.append(0)
            else:
                feature_list.append(1)

        return feature_list

    # 同时使用pf信息和压缩频谱信息
    def feature_creat(self, all_feature, linelen):
        feature = []
        strategy = 'full'
        if strategy == 'full':
            for key, value in all_feature.items():
                if key == "spectrum":
                    normalized_spectrum = [0 for _ in range(linelen)]
                    for line in value:
                        normalized_spectrum[line-1] = 1
                    feature += normalized_spectrum
                else:
                    feature += value

        return feature



# 特征距离的计算方法
class Distance:

    def usedfunction(self, v1, v2):
        # return self.kendall_tau_revised(v1, v2)
        # return self.get_single_mseer_distance(v1, v2)
        # return self.kendall_tau(v1, v2)
        return self.euclidean_space_distance(v1, v2)


    # 欧式空间距离
    def euclidean_space_distance(self, v1, v2):
        if not len(v1) == len(v2):
            print('空间向量长度不一致')
            return False
        # 距离
        D = 0
        # 空间维度
        dimensionality = len(v1)

        for i in range(dimensionality):
            D += math.pow((v1[i] - v2[i]), 2)
        D = math.pow(D, 0.5)
        return D

    def kendall_tau(self, v1, v2):
        # vi [sus1, sus2,...,susn]
        # 对多个语句进行的从大到小排序序列
        if not len(v1) == len(v2):
            print('序列向量不符')
            return False

        D = 0
        for i in range(len(v1)):
            for j in range(len(v2)):
                if i >= j:
                    continue
                if (v1[i]-v1[j])*(v2[i]-v2[j]) < 0:
                    D += 1

        return D

    def kendall_tau_revised(self, v1, v2):
        # vi [sus1, sus2,...,susn]
        # 对多个语句进行的从大到小排序序列
        if not len(v1) == len(v2):
            print('序列向量不符')
            return False

        D = 0
        for i in range(len(v1)):
            for j in range(len(v2)):
                if i >= j:
                    continue
                if (v1[i]-v1[j])*(v2[i]-v2[j]) < 0:
                    D += 1/v1[i] + 1/v1[j] + 1/v2[i] + 1/v2[j]

        return D

    def get_single_mseer_distance(self, list1, list2):
        code_num = len(list1)
        temp = 0
        for i in range(code_num):
            for j in range(code_num):
                if j <= i:
                    continue
                A = list1[i] - list1[j]
                b = list2[i] - list2[j]
                if A * b < 0:
                    # temp+=1
                    temp += 1 / list1[i] + 1 / list1[j] \
                            + 1 / list2[i] + 1 / list2[j]
                else:
                    continue
        return temp







class Kmeans:
    def __init__(self, data, clusternum=2):
        self.data = data
        self.clusternum = clusternum

    def ClusterMain(self):
        while True:
            if self.clusternum == 1:
                # print('无法聚类')
                # return False
                return {0: list(range(len(self.data)))}
            try:
                clasterfunction_out = sklearn_KMeans(n_clusters=self.clusternum).fit(self.data)
                break
            except:
                self.clusternum -= 1

        clusterdata = dict()
        for i in range(len(self.data)):
            label = clasterfunction_out.labels_[i]
            if label not in clusterdata:
                clusterdata[label] = []
            clusterdata[label].append(i)
        return clusterdata

class Agglomerative:
    def __init__(self, data, clusternum):
        self.data = data
        self.clusternum = clusternum

    def ClusterMain(self):
        try:
            clasterfunction_out = sklearn_AgglomerativeClustering(n_clusters=self.clusternum).fit(self.data)
        except:
            print('中心过多')
            return False

        clusterdata = dict()
        for i in range(len(self.data)):
            label = clasterfunction_out.labels_[i]
            if label not in clusterdata:
                clusterdata[label] = []
            clusterdata[label].append(i)
        return clusterdata


class MSeer:
    def __init__(self, data=[]):
        # print(self, data)
        self.file = ''
        self.alldata = data
        self.data = []
        self.cluster_sorteddata = dict()

        # 定义常量
        # self.theta = 7.93449    # 待定
        self.theta = 1.5    # 待定
        self.alpha = 4/(self.theta*self.theta)
        self.beta = 4/math.pow(1.5*self.theta, 2)
        # 潜在值阈值 大于等于0
        self.potential_threshold = 0
        self.M0 = -1


        # 【特征， 潜在值， 状态】
        self.potential_data = []
        # 中心
        self.center_loc_list = []

    def dataset(self):

        for point in self.alldata:
            if point not in self.data:
                self.data.append(point)
        return

    def settheta(self):
        # print(self.data)
        distance_list = []
        for i, point_i in enumerate(self.data):
            for j, point_j in enumerate(self.data):
                if j >= i:
                    continue
                distance_list.append(Distance().usedfunction(point_i, point_j))
        distance_list.sort()
        if len(distance_list) == 0:
            return
        cut = int(math.floor(len(distance_list) * 5 / 100))
        min_v = distance_list[cut]
        max_v = distance_list[len(distance_list) - 1 - cut]
        for i in range(cut):
            distance_list[i] = min_v
            distance_list[len(distance_list) - 1 - i] = max_v

        self.theta = sum(distance_list) / len(distance_list)/2
        return


    # MSeer
    def ClusterMain(self):
        # ---------重设数据---------
        # print(datetime.datetime.now(), 'start')
        self.dataset()
        # print(datetime.datetime.now(), 'dataset')
        self.settheta()
        # print(datetime.datetime.now(), 'settheta')
        # ---------设置中心点---------
        # 计算潜在值
        self.SetPotentialValue()
        # print(datetime.datetime.now(), '计算潜在值')

        # 循环计算中心
        while True:
            centerdata = self.FindNewCenter()
            if not centerdata:
                break
            self.ResetPotentialValue(centerdata)
        # print(datetime.datetime.now(), '找到中心')

        # ---------更新聚类及中心点---------
        medoidsnum = len(self.center_loc_list)
        sum_list = self.CenterSum(self.center_loc_list)

        initeration = True
        repeattime = 0
        while initeration:
            repeattime += 1
            if repeattime >= len(self.data)*2:
                break
            # print(datetime.datetime.now(), 'while', repeattime)
            initeration = False
            for i, [center_sum, cluster_data_list] in enumerate(sum_list):
                min_initeration_center_sum = sys.maxsize
                min_initeration_center_loc = -1

                for virtual_data_loc in cluster_data_list:
                    initeration_center_sum = self.IniterationCenterSum(virtual_data_loc, cluster_data_list)
                    # print('%sth center:%s sum:%s,virtualcenter:%s sum%s' % (i, self.center_loc_list[i], center_sum, virtual_data_loc, initeration_center_sum))
                    if initeration_center_sum < min_initeration_center_sum:
                        min_initeration_center_sum = initeration_center_sum
                        min_initeration_center_loc = virtual_data_loc
                # print(i, min_initeration_center_sum, center_sum)
                if min_initeration_center_sum < center_sum:
                    self.center_loc_list[i] = min_initeration_center_loc
                    # print(self.center_loc_list)
                    sum_list = self.CenterSum(self.center_loc_list)
                    initeration = True
                    break

        clusterdata = dict()
        for i, point in enumerate(self.alldata):
            for j, [sum, data] in enumerate(sum_list):
                if self.data.index(point) in data:
                    if j not in clusterdata:
                        clusterdata[j] = []
                    clusterdata[j].append(i)
        # print('1')
        self.cluster_sorteddata = dict()
        for key, value in clusterdata.items():
            D_metrix = []
            for fom1 in value:
                sumD = 0
                for fom2 in value:
                    if fom1 == fom2:
                        continue
                    sumD += Distance().usedfunction(self.alldata[fom1], self.alldata[fom2])
                D_metrix.append([fom1, sumD])
            # A = sorted(D_metrix, key=lambda x: x[1],)
            # B = list(map(lambda x: x[0], sorted(D_metrix, key=lambda x: x[1],)))
            self.cluster_sorteddata[key] = list(map(lambda x: x[0], sorted(D_metrix, key=lambda x: x[1],)))
        
        # print(f"聚类中心:{self.center_loc_list} len {len(self.center_loc_list)}")    
        return self.cluster_sorteddata





    # 初始计算潜在值后的data
    def SetPotentialValue(self):
        for i, point_i in enumerate(self.data):
            P = 0
            for j, point_j in enumerate(self.data):
                # print(i, math.pow(Distance().usedfunction(point_i, point_j), 2))
                P += math.pow(math.e, -1*self.alpha*math.pow(Distance().usedfunction(point_i, point_j), 2))
            self.potential_data.append([point_i, P, 'virtual'])
        return

    def find(self):
        r_theta = 7.93449
        r_alpha = 4/(r_theta*r_theta)
        potential_data = []
        for i, point_i in enumerate(self.alldata):
            P = 0
            for j, point_j in enumerate(self.alldata):
                # print(i, math.pow(Distance().usedfunction(point_i, point_j), 2))
                P += math.pow(math.e, -1*r_alpha*math.pow(Distance().usedfunction(point_i, point_j), 2))
            potential_data.append([point_i, P, 'virtual'])
        # print(r_theta, potential_data)
        # while r_theta < 7.94:
        #     r_theta += 0.00001
        #     r_alpha = 4/(r_theta*r_theta)
        #     potential_data = []
        #     for i, point_i in enumerate(self.alldata):
        #         P = 0
        #         for j, point_j in enumerate(self.alldata):
        #             # print(i, math.pow(Distance().usedfunction(point_i, point_j), 2))
        #             P += math.pow(math.e, -1*r_alpha*math.pow(Distance().usedfunction(point_i, point_j), 2))
        #         potential_data.append([point_i, P, 'virtual'])
        #     print(r_theta, potential_data)

    # 计算中心位置
    # 返回 中心loc 中心的特征值 中心的潜在值 不存在则返回 False
    def FindNewCenter(self):
        # 选定初始点
        maxloc_list = []
        maxvalue = 0
        for i, [point_i, potential, type] in enumerate(self.potential_data):
            if not type == 'virtual':
                continue
            if len(self.center_loc_list) > 0:
                if potential <= self.potential_threshold:
                    continue
            if potential > maxvalue:
                maxvalue = potential
                maxloc_list = []
                maxloc_list.append(i)
            elif potential == maxvalue:
                maxloc_list.append(i)

        while len(maxloc_list) > 0:
            rt_loc = maxloc_list.pop(random.sample(range(len(maxloc_list)), 1)[0])
            rt_R = self.potential_data[rt_loc][0]
            rt_M = self.potential_data[rt_loc][1]
            # 存在符合大于阈值的点
            if self.M0 == -1:
                # 选择M0
                self.M0 = rt_M
                return rt_loc, rt_R, rt_M
            else:
                # 选择M theta
                if rt_M > self.M0 * 0.5:
                    return rt_loc, rt_R, rt_M
                elif rt_M < self.M0 * 0.15:
                    return False
                else:
                    Dmin = sys.maxsize
                    for loc in self.center_loc_list:
                        distance = Distance().usedfunction(rt_R, self.data[loc])
                        if distance < Dmin:
                            Dmin = distance
                    if Dmin/self.theta + rt_M / self.M0 >= 1:
                        return rt_loc, rt_R, rt_M
                    else:
                        self.potential_data[rt_loc][1] = 0
        else:
            # 不存在大于阈值的点
            return False


    # 迭代更新潜在值
    def ResetPotentialValue(self, centerdata):
        self.potential_data[centerdata[0]][2] = 'center'
        self.center_loc_list.append(centerdata[0])
        for i, [point_i, potential, type] in enumerate(self.potential_data):
            if not type == 'virtual':
                continue
            newpotential = potential - centerdata[2]*\
                           math.pow(math.e, -1*self.beta*math.pow(Distance().usedfunction(point_i, centerdata[1]), 2))
            self.potential_data[i][1] = newpotential
        return

    # 给定中心点centerlist计算每个类簇的sum和其中的点
    def CenterSum(self, centerlist):
        sum_list = [[0, []] for _ in range(len(centerlist))]
        for i, virtualdata in enumerate(self.data):
            minsum = sys.maxsize
            minsum_loc = -1
            if i in centerlist:
                sum_list[centerlist.index(i)][1].append(i)
                continue
            for j, center_loc in enumerate(centerlist):
                distance = Distance().usedfunction(virtualdata, self.data[center_loc])
                if distance < minsum:
                    minsum = distance
                    minsum_loc = j
            sum_list[minsum_loc][0] += minsum
            sum_list[minsum_loc][1].append(i)
        return sum_list

    # 计算迭代中心和
    def IniterationCenterSum(self, center_loc, cluster_list):
        sum = 0
        for cluster_data_loc in cluster_list:
            sum += Distance().usedfunction(self.data[cluster_data_loc], self.data[center_loc])

        return sum

def read_json_file(file_path):
    """
    读取JSON文件并返回解析后的数据

    Parameters:
    - file_path (str): JSON文件的路径

    Returns:
    - data (dict): 解析后的JSON数据，如果文件不存在则返回None
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
        return None
    except Exception as e:
        print(f"发生错误: {e}")
        return None


def init(pid, vidList, pocessname):
    result = {}
    for vid in vidList:
      killJuzhenPath = f"/home/changzexing/MaxMajorKillJuzhen/{pid}/{vid}b/killJuzhen.json"
      if not os.path.exists(killJuzhenPath):
          print(f"{pid}-{vid} 不存在")
          continue
      # 第一步获取杀死矩阵
      killJuzhen = read_json_file(killJuzhenPath)
      if killJuzhen is None:
          continue
      datas = []
      mutantInfo = []
      for item in killJuzhen:
        datas.append(killJuzhen[item])
        mutantInfo.append(item)
    #   print(len(datas))
      myfile = f"/home/changzexing/MajorJuleiInfo/{pid}/{vid}b"
      if not os.path.exists(myfile):
        os.makedirs(myfile)
      with open(f"{myfile}/mutantInfo.json", 'w') as f:
        json.dump(mutantInfo, f)
      info = MSeer(datas).ClusterMain()
    #   result[f"{pid}-{vid}"] = info
    #   print(f"{pid}-{vid} 完成")
    #   print(MSeer(datas).ClusterMain())
    #   print(len(info), sum)
        # myfile = f"/home/changzexing/MajorJuleiInfo/{pid}"
        # if not os.path.exists(myfile):
        #     os.makedirs(myfile)
      with open(f"{myfile}/juleiInfo.json", 'w') as f:
        json.dump(info, f)
      print(f"{pid}完成")
    
def hunheinit(pid, vidList, pocessname):
    result = {}
    num = 0
    for vid in vidList:
        print(f'{pid} - {vid} 开始')
        majorKillJuzhenPath = f"/home/changzexing/MaxMajorKillJuzhen/{pid}/{vid}b/killJuzhen.json"
        mbertKillJuzhenPath = f"/home/changzexing/MaxMbertKillJuzhen/{pid}/{vid}b/killJuzhen.json"
        if not os.path.exists(majorKillJuzhenPath) and not os.path.exists(mbertKillJuzhenPath):
            print(f"{pid}-{vid} 杀死矩阵都不存在")
            continue
        # 第一步获取杀死矩阵
        majorKillJuzhenData = read_json_file(majorKillJuzhenPath)
        mbertKillJuzhenData = read_json_file(mbertKillJuzhenPath)
        if majorKillJuzhenData is None and mbertKillJuzhenData is None:
            print(f"{pid}-{vid} 杀死矩阵都为空")
            continue
        AllMajorMutant = {}
        AllMbertMutant = {}
        # print(len(majorKillJuzhenData))
        # print(len(mbertKillJuzhenData))
        if majorKillJuzhenData is not None:
            for fom in majorKillJuzhenData:
                fileLine = '-'.join(fom.split('-')[:-1])
                id = fom.split('-')[-1]
                if fileLine not in AllMajorMutant:
                    AllMajorMutant[fileLine] = [fom]
                else:
                    AllMajorMutant[fileLine].append(fom)
        if mbertKillJuzhenData is not None:
            for fom in mbertKillJuzhenData:
                fileLine = '-'.join(fom.split('-')[:-1])
                id = fom.split('-')[-1]
                if fileLine not in AllMbertMutant:
                    AllMbertMutant[fileLine] = [fom]
                else:
                    AllMbertMutant[fileLine].append(fom)
        # print(AllMajorMutant)
        # print(AllMbertMutant)
        # break
        # 获取两个变异工具每一行怀疑度最高的那个
        mutantResultPath = f"/home/changzexing/MajorMbertStatementSusType/{pid}/Ochiai.json"
        mutantResultData = read_json_file(mutantResultPath) 
        # mutantResultData格式：
        # "Chart-1b": {
        #   "source-org-jfree-chart-renderer-category-AbstractCategoryItemRenderer-235": {
        #     "major": {
        #       "rank": 145.5,
        #       "faulty": false,
        #       "sus": "0",
        #       "muttantInfo": {
        #         "mut_operator": "ROR",
        #         "bianhao": 2
        #       }
        #     },
        #     "mbert": {
        #       "rank": 138.0,
        #       "faulty": false,
        #       "sus": "0",
        #       "muttantInfo": {
        #         "mut_operator": "BinaryOperatorMutator",
        #         "orig_token": "&&",
        #         "pred_token": "||",
        #         "bainhao": "2"
        #       }
        #     }
        #   },
        # }
        pidMutantResult = mutantResultData[f"{pid}-{vid}b"]
        

        datas = []
        mutantInfo = []
        for item in pidMutantResult:
            # item: source-org-jfree-chart-renderer-category-AbstractCategoryItemRenderer-235
            # print(item)
            # print(pidMutantResult[item])
            if "major" in pidMutantResult[item] and "mbert" in pidMutantResult[item]:
                majorSus = pidMutantResult[item]["major"]["sus"]
                mbertSus = pidMutantResult[item]["mbert"]["sus"]
                # majorMutantName = f'{item}-{pidMutantResult[item]["major"]["muttantInfo"]["bianhao"]}'
                # mbertMutantName = f'{item}-{pidMutantResult[item]["mbert"]["muttantInfo"]["bainhao"]}'
                # majorMutantKillJuzhen = majorKillJuzhenData[majorMutantName]
                # mbertMutantKillJuzhen = mbertKillJuzhenData[mbertMutantName]
                if float(majorSus) > float(mbertSus):
                    if item in AllMajorMutant:
                        for tmp in AllMajorMutant[item]:
                            # datas.append(majorKillJuzhenData[tmp])
                            mutantInfo.append(f'major~{tmp}')
                else:
                    if item in AllMbertMutant:
                        for tmp in AllMbertMutant[item]:
                            # datas.append(mbertKillJuzhenData[tmp])
                            mutantInfo.append(f'mbert~{tmp}')
            elif "major" in pidMutantResult[item]:
                # majorMutantName = f'{item}-{pidMutantResult[item]["major"]["muttantInfo"]["bianhao"]}'
                # majorMutantKillJuzhen = majorKillJuzhenData[majorMutantName]
                # datas.append(majorMutantKillJuzhen)
                # mutantInfo.append(f'major~{majorMutantName}')
                if item in AllMajorMutant:
                    for tmp in AllMajorMutant[item]:
                        # datas.append(majorKillJuzhenData[tmp])
                        mutantInfo.append(f'major~{tmp}')
            elif "mbert" in pidMutantResult[item]:
                # mbertMutantName = f'{item}-{pidMutantResult[item]["mbert"]["muttantInfo"]["bainhao"]}'
                # mbertMutantKillJuzhen = mbertKillJuzhenData[mbertMutantName]
                # datas.append(mbertMutantKillJuzhen)
                # mutantInfo.append(f'mbert~{mbertMutantName}')
                if item in AllMbertMutant:
                    for tmp in AllMbertMutant[item]:
                        # datas.append(mbertKillJuzhenData[tmp])
                        mutantInfo.append(f'mbert~{tmp}')
        # print(len(mutantInfo))
        # if len(mutantInfo) == 0: num += 1
        # continue
        # break
        myfile = f"/home/changzexing/MajorMbertJuleiInfo/{pid}/{vid}b"
        if not os.path.exists(myfile):
            os.makedirs(myfile)
        with open(f"{myfile}/mutantInfo.json", 'w') as f:
            json.dump(mutantInfo, f)
        # info = MSeer(datas).ClusterMain()
        # with open(f"{myfile}/juleiInfo.json", 'w') as f:
        #     json.dump(info, f)
        print(f"{pid}-{vid} 完成")
        # print(num)
        # return  
        # myfile = f"/home/changzexing/MajorMbertJuleiInfo/{pid}"
        # if not os.path.exists(myfile):
        #     os.makedirs(myfile)
        
    print(f"{pid}完成")    

def startProcess( pid, svid, evid):
    num = 8
    No = []
    if pid == 'Chart':
       No = ['4', '5']
       num = evid - svid - len(No) + 1
    elif pid == 'Lang':
      No = []
    elif pid == 'Math':
      No = []
      num = 53
    elif pid == 'Time':
      No = ['21', '27']
      num = evid - svid - len(No) + 1
    elif pid == 'Mickito':
        No = ['18', '19', '20', '21', '35']
    len1 = math.ceil((evid - svid + 1 - int(len(No))) / num)
    # print(len1, len(No))
    versionList = {}
    for i in range(num):
      versionList[str(i)] = []
    tmp = 1
    no = 0
    for i in range(svid, evid + 1):
      if str(i) in No:
            continue
      if tmp == len1 + 1 and no < num - 1:
         tmp = 1
         no += 1
      tmp += 1
      versionList[str(no)].append(i)
    # print(versionList)
    # return
    # Create a process pool with 6 workers
    with ProcessPoolExecutor(max_workers=num) as pool:
        # Submit the jobs to the pool
        futures = [pool.submit(hunheinit, pid, versionList[str(i)], 'process%d'%(i,)) for i in range(num)]
        # futures = [pool.submit(init, pid, versionList[str(i)], 'process%d'%(i,)) for i in range(num)]
        
        # Wait for all the jobs to complete
        for future in futures:
            future.result()

if __name__ == '__main__':
    # init 处理major的聚类
    init("Chart", 1, 26)
    # init("Chart", range(1, 26), 'processname')
    # init("Closure", range(1, 5), 'processname')
    # init('Time', range(1, 26), 'processname')
    # init('Lang', range(1, 1), 'processname')
    # init('Mockito', range(1, 38), 'processname')
    # init('Math', 1, 106)
    # hunheinit 处理major和μbert结合的聚类
    # hunheinit("Chart", 1, 26)
    # hunheinit('Time', 1, 27)
    # hunheinit('Lang', 1, 65)
    # hunheinit('Math', 1, 106)
    # startProcess('Mockito', 12, 38)
    # startProcess('Cli', 1, 39)
    # startProcess('Codec', 1, 18)
    # startProcess('Csv', 1, 16) 
    # startProcess('Gson', 1, 18)
    # startProcess('JacksonXml', 1, 6)
    










