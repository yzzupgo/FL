import os
import shutil
import subprocess
import sys
import time
import json

import psutil


def projectPath():
    project_path = os.path.dirname(__file__)
    return project_path

class Datasave():
    def __init__(self):
        # 直接使用数据
        self.compile = False  # 是否可编译
        self.doc = ''  # 版本
        self.path = ''  # 被测程序版本路径
        self.src_tests = ''  # 测试用例路径
        self.Fault_Record = []  # 真实故障
        self.testcase = []  # 全体测试用例路径信息

        self.cov = []  # 原始程序频谱信息
        self.true_out = []  # 真实执行结果
        self.or_out = []  # 原始执行结果
        self.or_list = []  # 原始执行结果pf信息
        self.linelen = 0  # 原程序代码行数

        self.fom_list = []  # 一阶变异体全体信息 [self.fom]
        self.out_dic = {}  # {line:[原始执行输出]}
        self.fomnum = 0  # fom数量

        self.hom_list_all = []  # hom全体 [self.hom_i]
        self.hom_list = []  # 限制数量后的hom [self.hom_i] 不带out信息
        self.hom_out_list = []  # hom全体 [self.hom_i] 带out信息
        self.hom_out_dic = {}  # 各行hom执行结果对比信息列表 用于输入 get_touple {line: [out_or]}
        self.homnum = 0  # hom数量
        self.homs = []

        self.hom_num = {}  # hom数量
        self.claster_group = {}  # fom信息在各聚类方法下的类别 { method： {group: fom_i_message}}
        self.claster_hom_all = {}  # 各方法下的hom全体 { method： [self.hom_i]}
        self.claster_hom = {}  # 各方法下的限制数量后的hom { method： [self.hom_i]}

        self.time_spend = {
            'pre_del': 0,  # 开始执行前的数据处理时间
            'del': 0,  # 开始执行后的总体时间
            'getsus': 0,  # 获取怀疑度的列表的时间
            'hom_del': 0,  # hon阶段全体处理时间
        }

        # 被调用数据
        self.sbfl = {
            'cov': [],
            'time': 0}
        self.mbfl = {
            'sorce': [],
            'fomnum': 0,
            'homs': [],
        }
        self.fom = {
            'message': [[10, '+', '-', 'a += 1', 'a -= 1', '2'], [10, '+', '-', 'a += 1', 'a -= 1', '2']],  # 变异信息
            'spectrum': [],  # 频谱信息
            'out_list': [],  # 执行结果pf信息
            'out_or': [],  # 执行结果
            'claster': {},  # 聚类信息
            'time': 0,
        }
        self.hom = {
            'message': [],  # 变异信息
            'spectrum': [],  # 频谱信息
            'out_list': [],  # 执行结果pf信息
            'out_or': [],  # 执行结果
            'time': 0,  # 执行时间
            }


# 文件操作
class File:
    # 删除执行信息
    def remove_relative(self, program_path, if_fom=False):
        try:
            if if_fom:
                self.re(program_path)
            time.sleep(0.5)
            self.re(program_path[:-1] + 'exe')
            self.re(program_path[:-1] + 'gcno')
            self.re(program_path + '.gcov')
            self.re(program_path[:-1] + 'gcda')
            self.re(program_path[:-1] + 'out')
            self.re(program_path[:-1] + 'o')
            return 0
        except:
            print("remove relative error", program_path)

    # 删除文件
    def re(self, path):
        try:
            if os.path.exists(path):
                if sys.platform == 'linux':
                    subprocess.Popen("rm " + path, shell=True)
                else:
                    subprocess.Popen("del " + path, shell=True)
                    # print('删除成功 %s'%(path))
                time.sleep(0.5)
        except:
            print('删除失败', path)
        return 0

    def clear(self, path):
        try:
            del_list = os.listdir(path)
            for f in del_list:
                file_path = os.path.join(path, f)
                if os.path.isfile(file_path):
                    os.remove(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            return 0
        except:
            print("clear error", path)

    # 读取文件内容
    def read_line(self, file_name):
        try:
            with open(file_name, 'r', encoding='utf-8') as f:
                text = f.read().splitlines()
            return text
        except:
            print("读取" + file_name + "出错")
            return

    # 读取文件内容
    def read(self, file_name):
        try:
            with open(file_name, 'r', encoding='utf-8') as f:
                text = f.read()
            return text
        except:
            print("读取" + file_name + "出错")

    # 比较执行结果
    def compare_res(self, user_res, correct_res):
        flag = True
        try:
            if not os.path.exists(user_res):
                return False
            user_res_text = self.read_line(user_res)
            correct_res_text = self.read_line(correct_res)

            if len(user_res_text) != len(correct_res_text):
                flag = False
            if flag:
                for i in range(len(correct_res_text)):
                    if user_res_text[i] != correct_res_text[i]:
                        flag = False
        except:
            flag = False
        return flag


# 命令行语句
class Command:
    def __init__(self, program_path):
        self.timeout = 1.5
        # self.debug = True
        self.debug = False
        self.Path = os.path.dirname(os.path.abspath(program_path))
        self.filepath = os.path.abspath(program_path)
        self.program = os.path.basename(program_path)
        self.out = self.filepath[:-1] + 'out'
        self.exe = self.filepath[:-1] + 'exe'

    # 编译
    def compile(self, program_name):
        # 生成 .exe .gcno
        if sys.platform == 'linux':
            compile_command = 'cd %s && timeout %s gcc -fprofile-arcs -ftest-coverage %s -o %s -lm' % (
                self.Path, self.timeout, self.filepath, self.filepath[:-1] + 'exe')
        else:
            compile_command = 'cd %s && gcc -fprofile-arcs -ftest-coverage -o %s %s' % (self.Path, program_name, self.filepath)
        if self.debug:
            print('compile', compile_command)
        return compile_command

    # 运行
    def run(self, testdata):
        # 生成 .gcda .out
        if sys.platform == 'linux':
            run = 'cd %s && timeout -s SIGKILL %s %s <%s> %s' % (
                self.Path, self.timeout, self.exe, os.path.abspath(testdata),
                self.filepath[:-1] + 'out')
        else:
            run = 'cd %s && %s <%s> %s' % (self.Path, self.exe, os.path.abspath(testdata), self.out)
        if self.debug:
            print('run', run)
        # os.system(run)
        return run

    # 报告
    def report(self):
        if sys.platform == 'linux':
            report = 'cd %s && timeout %s gcov %s' % (self.Path, self.timeout, self.exe)
        else:
            report = 'cd %s && gcov %s' % (self.Path, self.program)
        if self.debug:
            print('report', report)
        return report


# 程序操作，包括编译文件，生成覆盖信息
class Programwork:
    def __init__(self):
        self.cov_result = [0, 0]
        self.cov_all = []
        self.res = []
        self.compile = False

    #  编译运行初始文件
    def program_start(self, program_path, testdata_dirpath):
        print(program_path, '***********************************开始程序***********************************')
        self.cov = False
        self.cov_result = [0, 0]
        command = Command(program_path)
        self.cov_all = []
        # 被执行行
        self.res = []
        # 执行情况
        for testdata in os.listdir(testdata_dirpath):
            if '.out'in testdata:
                continue
            print(testdata, '---------------testdata---------------')
            File().remove_relative(program_path)
            if not self.program_compile(program_path):
                return self.cov, self.cov_result
            info = []
            testdata_in = os.path.join(testdata_dirpath, testdata)
            testdata_out = testdata_in[:-2] + 'out'

            try:
                os.system(command.run(testdata_in))
            except:
                print('命令行运行错误')
                return self.cov, self.cov_result
            # print('system run end')

            try:
                os.system(command.report())
                while not os.path.exists(program_path + '.gcov'):
                    time.sleep(0.3)
            except:
                print('命令行报告错误')
                return self.cov, self.cov_result
            # print('system report end')

            try:
                # 读取覆盖信息
                lines = File().read_line(program_path + '.gcov')
                for line in lines:
                    # print(line)
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
                if len(info) == 0:
                    # print('输出错误', program_path)
                    return self.cov, self.cov_result
            except:
                # print('输出读取错误', program_path)
                return self.cov, self.cov_result

            res_i = File().compare_res(str(program_path[:-1]) + 'out', testdata_out)
            self.res.append(res_i)
            self.cov_all.append(info)
        self.cov = True
        File().remove_relative(program_path)
        # print('运行成功 ', program_path)
        return self.res, self.cov_all

    # 编译文件
    def program_compile(self, program_path):
        self.compile = False
        try:
            program_name = os.path.basename(program_path)[:-2]
            command = Command(program_path)

            p = subprocess.Popen(command.compile(program_name), stdout=subprocess.PIPE, shell=True, cwd=command.Path)
            p.communicate()
            if p.returncode:
                print('编译错误')
                return self.compile
            while not os.path.exists(os.path.join(command.Path, program_name + '.exe')):
                time.sleep(0.3)
            # print('编译成功')
            self.compile = True
            return self.compile
        except:
            return self.compile


# 通过tag文件获取真实故障信息
def tag_read_trueerror(tagpath):
    true_error = {}
    with open(tagpath, 'r') as f:
        text_lines = f.readlines()
    true_error['error_num'] = int(text_lines[2])
    true_error['error_lines'] = []
    for i in range(3, 3 + true_error['error_num']):
        true_error['error_lines'].append(text_lines[i].split('<TAG>')[1])
    return true_error


# 读取存储格式的变异体
def read_fom_message(info):
    [num, path, line, before, after] = info.split('----')
    return num, path, int(line), before, after


# 评估指标：
class Evaluation:

    # 根据怀疑度列表[[line_i,score]] 和真实错误[line_1]给出rank列表
    def rank(self, suslist, realfalt, line_len):
        fault_rank = []    #[[line, best, average, worst]
        suslist_sorted = sorted(suslist, key=lambda x: x[1], reverse=True)

        for fault in realfalt:
            fault_sus = -1
            r = 0
            tie = 0
            for [line, score] in suslist_sorted:
                if type(score) == list:
                    fault_sus = []
                if line == fault:
                    fault_sus = score
                    break
            if fault_sus == -1 or fault_sus == []:
                fault_rank.append([fault, len(suslist)+1, (len(suslist)+1 + line_len)/2, line_len])
                continue

            for i, [line, score] in enumerate(suslist_sorted):
                # if score > fault_sus:
                if sorted([score, fault_sus]).index(score) == 1 and not score == fault_sus:
                    r += 1
                elif score == fault_sus:
                    tie += 1
                else:
                    break
            fault_rank.append([fault, r+1, r+(max(tie, 1)+1)/2, r+max(tie, 1)])
        return fault_rank


def get_rank_value(ranks, linelen):
    resdata = {
        'rank_best': [],
        'rank_average': [],
        'rank_worst': [],
        'exam_best': [],
        'exam_average': [],
        'exam_worst': [],
    }
    rank = [[], [], []]
    for fault in ranks:
        rank[0].append(fault[1])
        rank[1].append(fault[2])
        rank[2].append(fault[3])
    resdata['rank_best'] = rank[0]
    resdata['rank_average'] = rank[1]
    resdata['rank_worst'] = rank[2]

    exam = [[], [], []]
    for ith in range(len(ranks)):
        if linelen > 0:
            exam[0].append(rank[0][ith]/linelen)
            exam[1].append(rank[1][ith]/linelen)
            exam[2].append(rank[2][ith]/linelen)
        else:
            exam[0].append(rank[0][ith])
            exam[1].append(rank[1][ith])
            exam[2].append(rank[2][ith])

    resdata['exam_best'] = exam[0]
    resdata['exam_average'] = exam[1]
    resdata['exam_worst'] = exam[2]
    return resdata


def json_rules(path):
    try:
        with open(path) as f_obj:
            data_json = json.load(f_obj)
    except:
        print('读取异常 %s' % path)
        return False
    return data_json

if __name__ == '__main__':
    print(projectPath())
    print(__file__)