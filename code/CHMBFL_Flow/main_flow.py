# Claster Based High Order Mutatiuon Based Fault Location
import os
import sys
import re
import util
import time
import datetime
from subprocess import Popen, PIPE
import subprocess
import numpy as np
import json
import random
from concurrent.futures import ThreadPoolExecutor

from codeflaws_version_control import Qr_excel
import sbfl.command as sbfl
import mbfl.command as mbfl

# 全局信息
import CHMBFL
from data_codeflaws import Codeflaws
# from CHMBFL_Flow.Cluster_Fom import cluster_sameline


class Codeflaws_old:
    def __init__(self):
        self.compile = False  # 是否可编译
        self.path = ''  # 被测程序路径
        self.testcase = []  # 测试用例信息
        self.true_out = []  # 真实执行结果
        self.Fault_Record = []  # 真实故障
        self.or_list = []  # 原始执行信息

        self.fom_list = []  # 变异信息 指示变异体如何变化
        self.out_list = []  # 执行信息 指示变异体执行结果
        self.fom_time = []  # 执行变异体用时
        self.out_dic = {}  # {line:[pf]}
        self.fomnum = 0  # fom数量
        self.time_spend = {
            'pre_del': 0,  # 开始执行前的数据处理时间
            'del': 0,  # 开始执行后的总体时间
            'getsus': 0  # 获取怀疑度的列表的时间
        }

        self.homnum = 0
        self.homs = []

    # 执行命令行
    def mysystem(self, cmd):
        try:
            ps = subprocess.Popen(cmd, shell=True, stdout=PIPE, stderr=subprocess.STDOUT)
            ps.wait()
            out, err = ps.communicate()
        except Exception as e:
            print(e)
        # print('\n\n %s \n\n' % out.decode("utf8","ignore"))
        if 'In function' in out.decode("utf8", "ignore") and 'error:' in out.decode("utf8", "ignore"):
            # print('编译错误', cmd)
            # print('编译错误', out.decode("utf8", "ignore"))
            return False
        elif 'Killed' in out.decode("utf8", "ignore"):
            # print('执行超时', out.decode("utf8", "ignore"))
            return False
        return True


    # 单独执行所有测试用例获取覆盖矩阵
    # [True, out, cov, test_src]
    def single_cov(self, name, src_path, tests_path):
        # sbfl使用 单独运行每个测试用例获取覆盖信息

        test_src = []  # 测试用例路径
        out = []  # 执行结果
        cov = []  # 测试用例覆盖信息
        gcov_path = src_path + '.gcov'
        true_com = util.Command(src_path)

        for i, test in enumerate(os.listdir(tests_path)):
            try:
                # 编译
                if not self.mysystem(true_com.compile(name)):
                    return [False]
                testcase = os.path.join(tests_path, test)
                test_src.append(testcase)

                st = datetime.datetime.now()
                # 运行
                self.mysystem(true_com.run(testcase))
                out.append(util.File().read(true_com.out))
                et = datetime.datetime.now()

                if i == 0:
                    ttime = (et - st).microseconds
                else:
                    ttime += (et - st).microseconds

                # 报告
                self.mysystem(true_com.report())
                time.sleep(0.3)
                res = self.read_gcov(gcov_path)
                if res[0]:
                    cov.append(res[1])
                else:
                    cov.append(False)
            except Exception as e:
                print(src_path, testcase, e)
                return [False]

        return [True, out, cov, test_src, ttime]

    # 判断编译是否通过
    def compile_pass(self, name, src_path):
        true_com = util.Command(src_path)
        # 编译文件
        if not self.mysystem(true_com.compile(name)):
            print('编译错误', src_path)
            return False
        return True

    # cov = False 不获取覆盖信息    return [True, test_src, true_out]
    # cov = all 获取全部覆盖信息
    # cov = list 执行失败测试用例覆盖信息
    # cov = i test_path 第i个测试用例
    def single(self, name, src_path, tests_path, cov=False):
        test_src = []  # 测试用例路径
        true_out = []  # 执行结果
        gcov_path = src_path + '.gcov'
        true_com = util.Command(src_path)

        # 编译文件
        if not self.mysystem(true_com.compile(name)):
            print('编译错误', src_path)
            return [False]

        if not cov:
            # 执行程序获取输出结果
            for test in os.listdir(tests_path):
                testcase = os.path.join(tests_path, test)
                test_src.append(testcase)
                if not self.mysystem(true_com.run(testcase)):
                    return [False]
                read_out = util.File().read(true_com.out)
                if not read_out:
                    return [False]
                true_out.append(read_out)
                os.remove(true_com.out)
            return [True, test_src, true_out]
        # 同时获取输出和全覆盖
        for i, test in enumerate(os.listdir(tests_path)):
            if type(cov) == list and cov[i] == 1:
                continue
            if type(cov) == int and not i == cov:
                continue
            testcase = os.path.join(tests_path, test)
            test_src.append(testcase)
            if not self.mysystem(true_com.run(testcase)):
                return [False]
            read_out = util.File().read(true_com.out)
            if not read_out:
                return [False]
            true_out.append(read_out)
            os.remove(true_com.out)

        self.mysystem(true_com.report())
        time.sleep(0.3)
        return_cov = self.read_gcov(gcov_path)
        return_cov.append(true_out)
        return return_cov

    # 读取gcov文件
    def read_gcov(self, path):
        info = []
        try:
            # 读取覆盖信息
            lines = util.File().read_line(path)
            for line in lines:
                line_cov = line.split(':', 2)
                line_cov[0] = line_cov[0].replace(' ', '')
                line_cov[1] = int(line_cov[1].replace(' ', ''))
                if line_cov[1] == 0:
                    continue
                try:
                    if int(line_cov[0]):
                        info.append(line_cov[1])
                except:
                    continue
                # print(line_cov)
            if len(info) == 0:
                print('输出错误', path)
                return [False, info]
        except:
            print('输出读取错误', path)
            return [False, info]
        return [True, info]

    # codeflaws 通过真实执行结果和原始执行结果对比出01向量
    def get_list_from_out(self, trueout, nowout):
        nowlist = []
        for i, out in enumerate(trueout):
            try:
                if len(out) == len(nowout[i]):
                    if out == nowout[i]:
                        nowlist.append(1)
                    else:
                        nowlist.append(0)
                else:
                    nowlist.append(0)
            except:
                nowlist.append(0)
        return nowlist


def sbfl_exec(src):
    try:
        print('执行sbfl')

        src_true = os.path.join(src, 'true_root', 'source', 'tcas.c')
        src_or = os.path.join(src, 'defect_root', 'source', 'tcas.c')
        src_tests = os.path.join(src, 'inputs')

        true_out = Codeflaws().single('tcas', src_true, src_tests)
        if not true_out[0]:
            print('获取真实信息失败')
            return [False]
        or_out = Codeflaws().single_cov('tcas', src_or, src_tests)
        if not or_out[0]:
            print('获取原始信息失败')
            return [False]
        outlist = Codeflaws().get_list_from_out(true_out[2], or_out[1])

        sbfl_touple = sbfl.GetTouleList(outlist, or_out[2])

        suslist = {}
        for line in sbfl_touple.keys():
            suslist[line] = sbfl_fun[0](sbfl_touple[line])

        sorted_sus = sorted(suslist.items(), key=lambda x: x[1], reverse=True)

        return [True, sorted_sus, or_out[4]]
    except Exception as e:
        print(e)
        return [False]


def single_fom_part(src, doc, Fault_Record, describe):
    # ----------------------------------------------------------------------------
    # 初始化数据
    data_return = util.Datasave()

    src_doc = os.path.join(src, doc, 'test_data')  # 题目大版本路径
    src_true = os.path.join(src_doc, 'true_root', 'source', 'tcas.c')  # 正确程序路径
    src_or = os.path.join(src_doc, 'defect_root', 'source', 'tcas.c')  # 故障程序路径
    src_tests = os.path.join(src_doc, 'inputs')  # 测试用例路径

    data_return.doc = doc
    data_return.path = src_or
    data_return.src_tests = src_tests
    data_return.Fault_Record = Fault_Record
    pre_del_s = datetime.datetime.now()

    # ----------------------------------------------------------------------------
    # 执行真实程序获取真实执行结果
    print('%s - 获取真实执行输出' % describe)
    singleout = Codeflaws().single('tcas', src_true, src_tests)
    if not singleout[0]:
        print('%s - 获取真实执行输出失败' % doc)
        return data_return
    else:
        # 全体测试用例信息， 正确程序执行输出
        data_return.testcase, data_return.true_out = singleout[1], singleout[2]

    # ----------------------------------------------------------------------------
    # 执行原始程序获取执行结果
    # print('%s - 获取原始执行输出' % describe)
    print('%s - 获取原始执行输出及频谱信息' % describe)
    singleout = Codeflaws().single_cov('tcas', src_or, src_tests)
    if not singleout[0]:
        print('%s - 获取原始执行输出及频谱信息失败' % describe)
        return data_return
    else:
        # 原始错误程序输出频谱
        # data_return.or_out = singleout[2]
        # data_return.or_list = Codeflaws().get_list_from_out(data_return.true_out, singleout[2])
        data_return.or_out = singleout[1]
        data_return.sbfl['cov'] = singleout[2]
        data_return.sbfl['time'] = singleout[4]
        data_return.or_list = Codeflaws().get_list_from_out(data_return.true_out, singleout[1])

    # ----------------------------------------------------------------------------
    # 获取执行行
    print('%s - 获取执行' % describe)
    singleout = Codeflaws().single('tcas', src_or, src_tests, data_return.or_list)
    if not singleout[0]:
        print('%s - 获取执行失败' % describe)
        # print('生成变异信息失败')
        return data_return
    # print(singleout)

    # ----------------------------------------------------------------------------
    # 生成一阶变异体
    mut_text = mbfl.Fom().fom_data('', src_or, singleout[1])
    muts = list(map(eval, mut_text.strip().split('\n')))

    # ----------------------------------------------------------------------------
    # 使用一阶变异数据进行全随机变异
    # hom_path = r'./report/c/Hom-%s-original.txt' % doc
    # hom_path = mbfl.Mutation().random_hmbfl(fompath, hom_path)

    # pre_del用时
    data_return.time_spend['pre_del'] = (datetime.datetime.now() - pre_del_s).microseconds
    # del用时开始
    del_s = datetime.datetime.now()

    # ----------------------------------------------------------------------------
    # 开始变异
    total_fom = len(muts)
    print('%s - 开始变异 \n fom num:' % describe, total_fom)

    for i, mut in enumerate(muts):
        if i % 10 == 0:
            print('%s -fom %s/%s/%s' % (datetime.datetime.now(), describe, i, total_fom))

        # 生成变异体文件
        src_mut = os.path.join(os.getcwd(), 'report', 'CHMBFL', 'fom', '%s-%s.c' % (doc, i))
        if not mbfl.Fom().mutation_build(src_or, src_mut, mut):
            continue

        # 执行变异体
        st = datetime.datetime.now()
        singleout = Codeflaws().single("%s-%s" % (doc, i), src_mut, src_tests, True)
        if not singleout[0]:
            continue
        et = datetime.datetime.now()

        # 获取输出信息
        out = Codeflaws().get_list_from_out(data_return.true_out, singleout[2])
        kill = Codeflaws().get_list_from_out(data_return.or_out, singleout[2])

        # **********保存数据
        fom = util.Datasave().fom
        fom['message'] = mut
        fom['spectrum'] = singleout[1]
        fom['out_list'] = out
        fom['kill_list'] = kill
        # fom['out_or'] = singleout[2]
        fom['time'] = (et - st).microseconds
        data_return.fom_list.append(fom)

    # del用时结束
    data_return.time_spend['del'] = (datetime.datetime.now() - del_s).microseconds

    data_return.fomnum = len(data_return.fom_list)

    # 可编译
    data_return.compile = True
    print('%s 执行一阶变异完成' % describe)
    return data_return


def single_hom_part(date_json, describe):
    # 使用的hom列表在文件中对应的key
    # hom_key = "hom_list"
    hom_key = "hom_list_all"

    # 数据预处理
    doc = list(date_json.keys())[0]
    data = date_json[doc]
    src_path = data['path']
    src_tests = data['src_tests']

    hom_list = data[hom_key]

    new_hom_list = []
    hom_out_list = []

    # ----------------------------------------------------------------------------
    # 开始变异
    del_s = datetime.datetime.now()
    total_hom = len(hom_list)

    print('%s - 开始变异 \n hom num:' % describe, total_hom)
    for i, mut in enumerate(hom_list):
        if i % 20 == 0:
            print('%s -hom %s/%s/%s' % (datetime.datetime.now(), describe, i, total_hom))

        # 跳过无法变异hom
        src_mut = os.path.join(os.getcwd(), 'report', 'CHMBFL', 'fom', '%s-%s.c' % (doc, i))
        if not mbfl.Fom().mutation_build(src_path, src_mut, mut['message']):
            continue

        # 执行变异体
        st = datetime.datetime.now()
        singleout = Codeflaws().single("%s-%s" % (doc, i), src_mut, src_tests, True)
        if not singleout[0]:
            continue
        et = datetime.datetime.now()

        # 获取输出信息
        out = Codeflaws().get_list_from_out(data['true_out'], singleout[2])
        kill = Codeflaws().get_list_from_out(data['or_out'], singleout[2])


        # **********保存数据
        new_hom_list.append(mut)
        mut['out_list'] = out
        # mut['out_or'] = singleout[2]
        mut['kill_list'] = kill
        mut['time'] = (et - st).microseconds
        hom_out_list.append(mut)

    data["time_spend"]["hom_del"] = (datetime.datetime.now() - del_s).microseconds
    data[hom_key] = new_hom_list
    data['hom_out_list'] = hom_out_list
    data['homnum'] = len(new_hom_list)
    date_json[doc] = data

    print('%s 执行二阶变异完成' % describe)
    return date_json


def hmbfl():
    # codeflaws
    src = os.path.join(os.getcwd(), 'cdata', 'version')
    print(src)

    completed = []
    if CHMBFL.continueread:
        for src in os.listdir(CHMBFL.doc_info_path):
            completed.append(src.split('_')[1])

    # 读取版本控制信息
    date_read = Qr_excel().read(CHMBFL.vc_path, CHMBFL.sheet)

    # 创建线程池
    executor = ThreadPoolExecutor(max_workers=CHMBFL.max_workers)
    task = []

    # 版本控制文件读取标记
    i = -1
    i_doc = {}

    while i < len(date_read) or len(task) > 0:
        if i >= 0:
            break
        while len(task) >= CHMBFL.max_workers:
            j = 0
            while j < len(task):
                if task[j].done():
                    t = task.pop(j)
                    data_return = t.result()
                    resdoc = data_return.doc
                    resi = i_doc[resdoc]
                    if not data_return.compile:
                        continue
                    # 获取真实故障
                    linelen = len(util.File().read_line(data_return.path))
                    # ----------------------------------------------------------------------------
                    # excel
                    date_read[resi, 3] = linelen
                    date_read[resi, 6] = data_return.fomnum
                    date_read[resi, 7] = data_return.homnum
                    date_read[resi, 15] = data_return.time_spend['pre_del']
                    date_read[resi, 16] = data_return.time_spend['del']

                    # ----------------------------------------------------------------------------
                    # json
                    date_json = {resdoc: {}}
                    date_json[resdoc]['sorce'] = data_return.or_list
                    date_json[resdoc]['fomnum'] = data_return.fomnum
                    date_json[resdoc]['homs'] = data_return.homs
                    date_json[resdoc]['homnum'] = data_return.homnum

                    date_json[resdoc]['realfault'] = data_return.Fault_Record
                    date_json[resdoc]['compile'] = data_return.compile
                    date_json[resdoc]['path'] = data_return.path
                    date_json[resdoc]['true_out'] = data_return.true_out
                    date_json[resdoc]['or_list'] = data_return.or_list
                    date_json[resdoc]['out_dic'] = data_return.out_dic
                    date_json[resdoc]['time_spend'] = data_return.time_spend
                    date_json[resdoc]['testcase'] = data_return.testcase
                    date_json[resdoc]['linelen'] = linelen

                    # ********存储数据
                    Qr_excel().save(CHMBFL.vc_path, CHMBFL.sheet, date_read)
                    json_path = os.path.join(CHMBFL.doc_info_path, CHMBFL.doc_info_name % resdoc)
                    with open(json_path, 'w') as f_obj:
                        json.dump(date_json, f_obj)

                else:
                    j += 1
                time.sleep(5)

        if i < len(date_read):
            i += 1
            if date_read[i, 2] == 0:
                # 跳过无法编译文件
                continue

            # 跳过continueread的数据
            doc = date_read[i, 0]
            Fault_Record = eval(date_read[i, 1])
            i_doc[doc] = i

            if doc in completed:
                continue
            describe = str(doc) + '-' + str(i)
            print(describe)
            # ----------------------------------------------------------------------------
            task.append(executor.submit(single_fom_part, src, doc, Fault_Record, describe))


# 生成一阶变异体
def testfom():
    # codeflaws  test

    date_json = {}

    src = os.path.join(os.getcwd(), 'testdata')
    print(src)

    doc = 'v1739'
    Fault_Record = [39, 40]
    fom_Data = single_fom_part(src, doc, Fault_Record, 'test')
    linelen = len(util.File().read_line(fom_Data.path))

    date_json[doc] = {}
    for name, value in vars(fom_Data).items():
        date_json[doc][name] = value
    date_json[doc]['linelen'] = linelen

    json_path = os.path.join(CHMBFL.doc_info_path, CHMBFL.doc_info_name % doc)
    with open(json_path, 'w') as f_obj:
        json.dump(date_json, f_obj)
    print('%s 一阶信息生成完成' % doc)


# 生成一阶变异体
def testhom():
    # codeflaws  test

    date_json = {}

    src = os.path.join(os.getcwd(), 'cdata', 'version')
    print(src)

    # --------------------------------------------------------------------------------
    # 数据初始化
    doc = 'v1739'
    Fault_Record = [39, 40]
    json_path = os.path.join(CHMBFL.doc_info_path, CHMBFL.doc_info_name % doc)

    # --------------------------------------------------------------------------------
    # 执行一阶变异体
    fom_Data = single_fom_part(src, doc, Fault_Record, 'test')
    if not fom_Data.compile:
        return False
    else:
        # 获取执行数据
        linelen = len(util.File().read_line(fom_Data.path))
        date_json[doc] = {}
        for name, value in vars(fom_Data).items():
            date_json[doc][name] = value
        date_json[doc]['linelen'] = linelen
        # 存储一阶执行获得的信息
        with open(json_path, 'w') as f_obj:
            json.dump(date_json, f_obj)
        print('%s 生成完成' % doc)

    # --------------------------------------------------------------------------------
    # 对一阶变异体进行聚类
    cluster_return = cluster(date_json)
    if not cluster_return:
        print('%s-聚类失败' % doc)
        return False
    else:
        date_json = cluster_return
        with open(json_path, 'w') as f_obj:
            json.dump(date_json, f_obj)
        print('%s-聚类完成' % doc)

    # --------------------------------------------------------------------------------
    # 使用聚类后的一阶变异体生成二阶变异体
    hom_return = mbfl.Mutation().claster_hom_random(date_json)
    if not hom_return:
        print('%s-二阶变异体生成失败' % doc)
        return False
    else:
        date_json = hom_return
        with open(json_path, 'w') as f_obj:
            json.dump(date_json, f_obj)
        print('%s-二阶变异体生成完成' % doc)

    # --------------------------------------------------------------------------------
    # 执行二阶变异体
    hom_out_return = single_hom_part(date_json, 'test')
    if len(hom_out_return[doc]["hom_out_list"]) == 0:
        print('%s-二阶变异体执行失败' % doc)
        return False
    else:
        date_json = hom_out_return
        with open(json_path, 'w') as f_obj:
            json.dump(date_json, f_obj)
        print('%s-二阶变异体执行完成' % doc)

    return


def single_hom(json_path, src, doc, Fault_Record, describe='undefine'):
    date_json = dict()
    # --------------------------------------------------------------------------------
    # 执行一阶变异体
    fom_Data = single_fom_part(src, doc, Fault_Record, describe)
    if not fom_Data.compile:
        return describe+'False'
    else:
        # 获取执行数据
        linelen = len(util.File().read_line(fom_Data.path))
        date_json[doc] = {}
        for name, value in vars(fom_Data).items():
            date_json[doc][name] = value
        date_json[doc]['linelen'] = linelen
        # 存储一阶执行获得的信息
        # with open(json_path, 'w') as f_obj:
        #     json.dump(date_json, f_obj)
        print('%s-一阶变异体生成完成' % doc)

    # --------------------------------------------------------------------------------
    # 生成二阶变异体全集
    # hom_return = mbfl.Mutation().all_hom_random(date_json)
    hom_return = mbfl.Mutation().all_nshom_random(date_json)
    if not hom_return:
        print('%s-二阶变异体全集生成失败' % doc)
        return describe+'False'
    else:
        date_json = hom_return
        # with open(json_path, 'w') as f_obj:
        #     json.dump(date_json, f_obj)
        print('%s-二阶变异体全集生成完成' % doc)

    # if len(hom_return[doc]['hom_list_all']) < 5000 or len(hom_return[doc]['hom_list_all']) > 16000:
    #     print('%s-二阶变异体数量不符' % doc)
    #     return describe+'False'


    # --------------------------------------------------------------------------------
    # 执行二阶变异体
    hom_out_return = single_hom_part(date_json, describe)
    if not hom_out_return:
        return False
    if len(hom_out_return[doc]["hom_out_list"]) == 0:
        print('%s-二阶变异体执行失败' % doc)
        return describe+'False'
    else:
        date_json = hom_out_return
        with open(json_path, 'w') as f_obj:
            json.dump(date_json, f_obj)
        print('%s-二阶变异体执行完成' % doc)
    return describe+'True'


def main():
    # codeflaws
    src = os.path.join(os.getcwd(), 'cdata', 'version')
    nsrc = os.path.join(os.getcwd(), 'Example')
    print(nsrc)

    # 读取版本控制信息
    # date_read = Qr_excel().read(CHMBFL.vc_path, CHMBFL.sheet)
    date_read = Qr_excel().read_old(CHMBFL.vc_path, CHMBFL.sheet)

    # 创建线程池
    executor = ThreadPoolExecutor(max_workers=CHMBFL.max_workers)
    task = []

    # 版本控制文件读取标记
    i = -1
    i_doc = {}
    # ip_list = ['202.4.157.11', '202.4.157.19', '222.199.230.227', '202.4.130.29']
    ip_list = ['202.4.157.11', '202.4.157.19', '222.199.230.148', '202.4.130.29']
    while i < len(date_read) - 1:
        i += 1
        if len(os.listdir(CHMBFL.doc_info_path)) >= 40:
            print('enough')
            break
        else:
            print('rest %s' % str(len(os.listdir(CHMBFL.doc_info_path))))


        # if not i % len(ip_list) == ip_list.index(CHMBFL.get_host_ip()):
        #     continue
        doc = date_read[i, 0]
        i_doc[doc] = i

        if CHMBFL.continueread:
            havedoc = False
            for file in os.listdir(CHMBFL.doc_info_path):
                if doc in file:
                    havedoc = True
            if havedoc:
                continue


        describe = str(doc) + '-' + str(i)
        # ----------------------------------------------------------------------------
        # 获取真实故障
        Fault_Record = eval(date_read[i, 1])
        json_path = os.path.join(CHMBFL.doc_info_path, CHMBFL.doc_info_name % doc)
        single_hom(json_path, src, doc, Fault_Record, describe)
        # print('放入%s' % doc)


def new_main():
    # codeflaws
    nsrc = os.path.join(os.getcwd(), 'Example')
    doc = 'v1647'
    Fault_Record = [11, 17]
    print(nsrc)

    # 读取版本控制信息
    # date_read = Qr_excel().read(CHMBFL.vc_path, CHMBFL.sheet)
    # date_read = Qr_excel().read_old(CHMBFL.vc_path, CHMBFL.sheet)

    # 创建线程池
    executor = ThreadPoolExecutor(max_workers=CHMBFL.max_workers)
    task = []

    # 版本控制文件读取标记

    # if len(os.listdir(CHMBFL.doc_info_path)) >= 40:
    #     print('enough')
    #     return
    # else:
    #     print('rest %s' % str(len(os.listdir(CHMBFL.doc_info_path))))

    describe = doc
    # ----------------------------------------------------------------------------
    # 获取真实故障
    json_path = os.path.join(r'/home/cyxy/access/Example', 'Fom_n%s_Feature.json' % doc)
    single_hom(json_path, nsrc, doc, Fault_Record, describe)
    # print('放入%s' % doc)
    return


def single_fom_part_1():
    describe = 'v1735'
    src = '../cdata/version'
    doc = 'v1750'
    # ----------------------------------------------------------------------------
    # 初始化数据
    data_return = util.Datasave()

    src_doc = os.path.join(src, doc, 'test_data')  # 题目大版本路径
    src_true = os.path.join(src_doc, 'true_root', 'source', 'tcas.c')  # 正确程序路径
    src_or = os.path.join(src_doc, 'defect_root', 'source', 'tcas.c')  # 故障程序路径
    src_tests = os.path.join(src_doc, 'inputs')  # 测试用例路径

    data_return.doc = doc
    data_return.path = src_or
    data_return.src_tests = src_tests
    pre_del_s = datetime.datetime.now()

    # ----------------------------------------------------------------------------
    # 执行真实程序获取真实执行结果
    print('%s - 获取真实执行输出' % describe)
    singleout = Codeflaws().single('tcas', src_true, src_tests)
    if not singleout[0]:
        print('%s - 获取真实执行输出失败' % doc)
        return data_return
    else:
        # 全体测试用例信息， 正确程序执行输出
        data_return.testcase, data_return.true_out = singleout[1], singleout[2]

    # ----------------------------------------------------------------------------
    # 执行原始程序获取执行结果
    # print('%s - 获取原始执行输出' % describe)
    # print('%s - 获取原始执行输出及频谱信息' % describe)
    # singleout = Codeflaws().single_cov('tcas', src_or, src_tests)
    # if not singleout[0]:
    #     print('%s - 获取原始执行输出及频谱信息失败' % describe)
    #     return data_return
    # else:
    #     # 原始错误程序输出频谱
    #     # data_return.or_out = singleout[2]
    #     # data_return.or_list = Codeflaws().get_list_from_out(data_return.true_out, singleout[2])
    #     data_return.sbfl['cov'] = singleout[2]
    #     data_return.sbfl['time'] = singleout[4]
    #     data_return.or_list = Codeflaws().get_list_from_out(data_return.true_out, singleout[1])
    data_return.or_out = ['5 1', '1 2', '7 8', '12 13', '16 17', '1 2', '1 2', '86 87', '7 8', '1 2', '81 82', '1 2', '36 37', '2 3', '1 2', '1 2', '3 4', '1 2', '2 3', '4 1', '1 2', '5 1', '1 2', '2 3', '1 2', '1 2', '3 4']
    data_return.true_out = ['5 1', '1 2', '7 8', '12 13', '16 17', '1 2', '1 2', '86 87', '7 8', '1 2', '81 82', '1 2', '36 37', '2 3', '1 2', '1 2', '3 4', '1 0', '2 0', '4 1', '1 0', '5 1', '1 2', '2 3', '1 2', '1 2', '3 4']

    # ----------------------------------------------------------------------------
    # 获取执行行
    print('%s - 获取执行' % describe)
    singleout = Codeflaws().single('tcas', src_or, src_tests, data_return.or_list)
    if not singleout[0]:
        print('%s - 获取执行失败' % describe)
        # print('生成变异信息失败')
        return data_return
    # print(singleout)

    # ----------------------------------------------------------------------------
    # 生成一阶变异体
    # mut_text = mbfl.Fom().fom_data('', src_or)
    # muts = list(map(eval, mut_text.strip().split('\n')))

    muts = [
        [[34, '    if( f > n ) f = f%n ;', '    if( f > 3 ) f = f%n ;', '?', '?', 10]],
        [[34, '    if( f > n ) f = f%n ;', '    if( f > 2 ) f = f%n ;', '?', '?', 10]],
        [[34, '    if( f > n ) f = f%n ;', '    if( f > 1 ) f = f%n ;', '?', '?', 10]],
    ]

    # ----------------------------------------------------------------------------
    # 使用一阶变异数据进行全随机变异
    # hom_path = r'./report/c/Hom-%s-original.txt' % doc
    # hom_path = mbfl.Mutation().random_hmbfl(fompath, hom_path)

    # pre_del用时
    data_return.time_spend['pre_del'] = (datetime.datetime.now() - pre_del_s).microseconds
    # del用时开始
    del_s = datetime.datetime.now()

    # ----------------------------------------------------------------------------
    # 开始变异
    total_fom = len(muts)
    print('%s - 开始变异 \n fom num:' % describe, total_fom)

    for i, mut in enumerate(muts):
        if i % 10 == 0:
            print('%s -fom %s/%s/%s' % (datetime.datetime.now(), describe, i, total_fom))

        # 生成变异体文件
        src_mut = os.path.join(os.getcwd(), '../report', 'CHMBFL', 'fom', '%s-%s.c' % (doc, i))
        if not mbfl.Fom().mutation_build(src_or, src_mut, mut):
            continue

        # 执行变异体
        st = datetime.datetime.now()
        singleout = Codeflaws().single("%s-%s" % (doc, i), src_mut, src_tests, True)
        if not singleout[0]:
            continue
        et = datetime.datetime.now()

        # 获取输出信息
        out = Codeflaws().get_list_from_out(data_return.true_out, singleout[2])
        kill = Codeflaws().get_list_from_out(data_return.or_out, singleout[2])

        # **********保存数据
        fom = util.Datasave().fom
        fom['message'] = mut
        fom['spectrum'] = singleout[1]
        fom['out_list'] = out
        fom['kill_list'] = kill
        # fom['out_or'] = singleout[2]
        fom['time'] = (et - st).microseconds
        data_return.fom_list.append(fom)

    # del用时结束
    data_return.time_spend['del'] = (datetime.datetime.now() - del_s).microseconds

    data_return.fomnum = len(data_return.fom_list)

    # 可编译
    data_return.compile = True
    print('%s 执行一阶变异完成' % describe)
    return data_return

if __name__ == '__main__':
    single_fom_part_1()
