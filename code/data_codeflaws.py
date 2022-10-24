import os
import sys
import re
import util
# import sbfl.command as sbfl_c
# import mbfl.command as mbfl_c
import time
import datetime
from subprocess import Popen, PIPE
import subprocess
import numpy as np
from codeflaws_version_control import Qr_excel
import json
import random



class Codeflaws:
    def __init__(self):
        self.compile = False  # 是否可编译
        self.path = ''  # 被测程序路径
        self.testcase = []  # 测试用例信息
        self.true_out = []  # 真实执行结果
        self.Fault_Record = []  # 真实故障
        self.or_list = []  # 原始执行信息

        self.fom_list = []  # 变异信息 指示变异体如何变化
        self.out_list = []  # 执行信息 指示变异体执行结果
        self.fom_time = [] #执行变异体用时
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
            return False
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
            print('%s test:%s-%s-%s' % (datetime.datetime.now(), test, i, len(os.listdir(tests_path))))
            util.File().remove_relative(src_path)
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
                    ttime = (et-st).microseconds
                else:
                    ttime += (et-st).microseconds

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
    def single_old(self, name, src_path, tests_path, cov=False):
        pout = False
        test_src = []  # 测试用例路径
        true_out = []  # 执行结果
        gcov_path = src_path + '.gcov'
        true_com = util.Command(src_path)

        # print('%s 开始编译' % datetime.datetime.now())
        # 编译文件
        if not self.mysystem(true_com.compile(name)):
            if pout:print('编译错误', src_path)
            return [False]

        if not cov:
            # 执行程序获取输出结果
            for test in os.listdir(tests_path):
                # print('%s 开始测试 %s' % (datetime.datetime.now(), test))

                testcase = os.path.join(tests_path, test)
                test_src.append(testcase)
                if not self.mysystem(true_com.run(testcase)):
                    if pout:print('运行错误 %s %s ' % (src_path, testcase))
                    return [False]
                read_out = util.File().read(true_com.out)
                if not read_out:
                    if pout:print('读取错误 %s' % true_com.out)
                    return [False]
                true_out.append(read_out)
                os.remove(true_com.out)
            return [True, test_src, true_out]

        for i, test in enumerate(os.listdir(tests_path)):
            if type(cov) == list and cov[i] == 1:
                continue
            if type(cov) == int and not i == cov:
                continue
            testcase = os.path.join(tests_path, test)
            self.mysystem(true_com.run(testcase))

        self.mysystem(true_com.report())
        time.sleep(0.3)
        return self.read_gcov(gcov_path)

    # cov = False 不获取覆盖信息    return [True, test_src, true_out]
    # cov = all 获取全部覆盖信息
    # cov = list 执行失败测试用例覆盖信息
    # cov = i test_path 第i个测试用例
    def single(self, name, src_path, tests_path, cov=False):
        pout = False
        test_src = []  # 测试用例路径
        true_out = []  # 执行结果
        gcov_path = src_path + '.gcov'
        true_com = util.Command(src_path)

        # 编译文件
        if not self.mysystem(true_com.compile(name)):
            if pout:print('编译错误', src_path)
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


    # cov = False 不获取覆盖信息    return [True, test_src, true_out]
    # cov = all 获取全部覆盖信息
    # cov = list 执行失败测试用例覆盖信息
    # cov = i test_path 第i个测试用例
    def single_testlist(self, name, src_path, tests):
        pout = False
        test_src = []  # 测试用例路径
        true_out = []  # 执行结果
        gcov_path = src_path + '.gcov'
        true_com = util.Command(src_path)

        # print('%s 开始编译' % datetime.datetime.now())
        # 编译文件
        if not self.mysystem(true_com.compile(name)):
            if pout:print('编译错误', src_path)
            return [False]

        # 执行程序获取输出结果
        for testcase in tests:
            # print('%s 开始测试 %s' % (datetime.datetime.now(), test))

            test_src.append(testcase)
            if not self.mysystem(true_com.run(testcase)):
                if pout:print('运行错误 %s %s ' % (src_path, testcase))
                return [False]
            read_out = util.File().read(true_com.out)
            if not read_out:
                if pout:print('读取错误 %s' % true_com.out)
                return [False]
            true_out.append(read_out)
            os.remove(true_com.out)
        return [True, test_src, true_out]



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
                # print('输出错误', path)
                return [False, info]
        except:
            # print('输出读取错误', path)
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

# 变异规则
def MutationRule(message1, meggage2):
        rule = 'Normal'
        # 常规变异
        if rule == 'Normal':
            # 不对同一行进行变异
            if message1[0] == meggage2[0]:
                return False
            else:
                return True
        # 非叠加变异
        if rule == 'NonSuperposed':
            # 不同行可变异
            # 同一行 其中一个是delet不能变异
            # 同一行算子不同可以变异
            if not message1[0] == meggage2[0]:
                return True
            elif message1[3] == 'delet' or meggage2[3] == 'delet':
                return False
            elif not message1[3] == meggage2[3] or \
                    (message1[5] < meggage2[5] and message1[5]+len(message1[3]) > meggage2[5]) or \
                    (meggage2[5] < message1[5] and meggage2[5]+len(meggage2[3]) > message1[5]):
                return True
            else:
                return False


