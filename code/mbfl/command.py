import datetime
import math
import os
import random

import util
import mbfl.mbfl_for as mbfl_for
from data_codeflaws import Codeflaws
from mbfl.mutpolyn import java_mutation
from mbfl.mutpolyn import mutation_trick
from mbfl.mutpolyn import new_mutation_trick
from util import File





# 变异体策略
class Mutation:
    # numset为每个类中包含的fom数量，n为需要的总数,为一个升序列表
    # 该算法求每个类中需要挑选的数量
    def numset(self, numlist, n):
        claster_num = len(numlist)
        ave = math.ceil((n * 2 / (claster_num * (claster_num - 1))) ** (1 / 2))
        if numlist[0] < ave:
            return self.numset(numlist[1:], math.ceil(n / numlist[0]))
        else:
            return [ave for i in range(claster_num)]

    # 新文件模式---------------------------------------------------------------------
    # 按照ntimeshom倍fom数量生成hom
    # 存在不同类簇数量差异大的问题
    def claster_hom_random(self, data_json, ntimeshom=5):
        doc = list(data_json.keys())[0]
        data = data_json[doc]
        if "claster_group" not in data:
            print("未进行聚类")
            return False

        claster_group = data["claster_group"]
        fomnum = data["fomnum"]
        print('%s 开始生成聚类二阶变异体' % datetime.datetime.now())

        # 相同类簇组合生成二阶变异体全体:
        claster_hom = {}
        for method in claster_group.keys():
            claster_hom[method] = []
            claster_method = claster_group[method].keys()
            for group in claster_method:
                for i, fom_i in enumerate(claster_group[method][group]):
                    for j, fom_j in enumerate(claster_group[method][group]):
                        if j <= i or fom_i[0][0] == fom_j[0][0]:
                            continue
                        hom = util.Datasave().hom
                        hom['message'].append(fom_i[0])
                        hom['message'].append(fom_j[0])
                        claster_hom[method].append(hom)
        data['claster_hom_all'] = claster_hom

        hom_list = []
        # 随机删除多余二阶变异体
        hom_num = {}
        for method in claster_hom.keys():
            while len(claster_hom[method]) > fomnum * ntimeshom:
                position = random.sample(range(len(claster_hom[method])), 1)
                claster_hom[method].pop(position[0])
            hom_num[method] = len(claster_hom[method])
            for hom in claster_hom[method]:
                if not hom in hom_list:
                    hom_list.append(hom)

        data['claster_hom'] = claster_hom
        data['hom_list'] = hom_list
        data['hom_out_list'] = []
        data['hom_num'] = hom_num
        data_json[doc] = data
        return data_json

    # 新文件模式---------------------------------------------------------------------
    # 变异体数量全生成
    def all_hom_random(self, data_json):
        doc = list(data_json.keys())[0]
        data = data_json[doc]

        fom_list = data["fom_list"]
        fomnum = data["fomnum"]
        print('%s 开始生成全二阶变异体' % datetime.datetime.now())

        # 全二阶变异体:
        hom_list_all = []
        for i, fom_i in enumerate(fom_list):
            for j, fom_j in enumerate(fom_list):
                if j <= i:
                    continue
                if fom_i['message'][0][0] == fom_j['message'][0][0] and fom_i['message'][0][3] == fom_j['message'][0][3] :
                    continue
                hom = util.Datasave().hom
                hom['message'].append(fom_i['message'][0])
                hom['message'].append(fom_j['message'][0])
                hom_list_all.append(hom)

        data['hom_list_all'] = hom_list_all
        data_json[doc] = data
        return data_json


    # 新文件模式---------------------------------------------------------------------
    # 变异体数量全生成
    def all_nshom_random(self, data_json):
        doc = list(data_json.keys())[0]
        data = data_json[doc]

        fom_list = data["fom_list"]
        fomnum = data["fomnum"]
        print('%s 开始生成全二阶变异体' % datetime.datetime.now())

        # 非叠加二阶变异体:
        hom_list_all = []
        for i, fom_i in enumerate(fom_list):
            for j, fom_j in enumerate(fom_list):
                if j <= i:
                    continue
                if not fom_i['message'][0][0] == fom_j['message'][0][0]:
                    continue
                hom = util.Datasave().hom
                hom['message'].append(fom_i['message'][0])
                hom['message'].append(fom_j['message'][0])
                hom_list_all.append(hom)

        data['hom_list_all'] = hom_list_all
        data_json[doc] = data
        return data_json


    # 旧文件模式---------------------------------------------------------------------
    # 使用筛选后的fom信息生成随机二阶变异体
    def random_hmbfl(self, fom_path, hom_path=''):
        if hom_path == '':
            hom_path = r'./report/c/Hom-original.txt'
        top = 1000
        fom_record = ''
        lines = util.File().read_line(fom_path)
        fom_dic = {}  # line: [[totalfom], [randomfom]]
        # fom按行分组
        for i, fommessage1 in enumerate(lines):
            fom1 = eval(fommessage1)[0]
            if fom1[0] not in fom_dic:
                fom_dic[fom1[0]] = [[], []]
            fom_dic[fom1[0]][0].append(fom1)
        # 每行机挑选fom
        for line in fom_dic:
            for i in random.sample(range(0, len(fom_dic[line][0])), min(len(fom_dic[line][0]), top)):
                fom_dic[line][1].append(fom_dic[line][0][i])

        for line1 in fom_dic:
            for line2 in fom_dic:
                if line1 >= line2:
                    continue
                for fom1 in fom_dic[line1][1]:
                    for fom2 in fom_dic[line2][1]:
                        if fom1 == fom2:
                            continue
                        fom_record += str([fom1, fom2]) + '\n'
        with open(os.path.join(hom_path), 'w', encoding='utf-8') as f:
            f.write(fom_record)
        return hom_path

    # 分集合组合策略    集合coll为行列表
    def Portfolio_hmbfl(self, fom_path, coll1, coll2):
        print('开始生成二阶变异体')
        hom_path = r'./report/c/Hom-original.txt'
        top = 1  # 每条语句使用的变异体数量
        fom_record = ''

        strategy = []  # 存储以生成的组合
        coll1_line = {}  # 存储coll1语句使用的变异体次数
        coll2_line = {}  # 存储coll2语句使用的变异体次数

        foms = util.File().read_line(fom_path)
        for i, fommessage1 in enumerate(foms):
            fom1 = eval(fommessage1)[0]
            if not fom1[0] in coll1:
                continue
            if not fom1[0] in coll1_line:
                coll1_line[fom1[0]] = []
            else:
                if len(coll1_line[fom1[0]]) >= top and i not in coll1_line[fom1[0]]:
                    continue

            for j, fommessage2 in enumerate(foms):
                fom2 = eval(fommessage2)[0]
                if not fom2[0] in coll2:
                    continue
                if not fom2[0] in coll2_line:
                    coll2_line[fom2[0]] = []
                else:
                    if len(coll2_line[fom2[0]]) >= top and j not in coll2_line[fom2[0]]:
                        continue

                if i not in coll1_line[fom1[0]]:
                    coll1_line[fom1[0]].append(i)
                if j not in coll2_line[fom2[0]]:
                    coll2_line[fom2[0]].append(j)
                if fom1[0] == fom2[0]:
                    continue
                if [i, j] in strategy or [j, i] in strategy:
                    continue
                strategy.append([i, j])
                fom_record += str([fom1, fom2]) + '\n'

        with open(os.path.join(hom_path), 'w', encoding='utf-8') as f:
            f.write(fom_record)
        print('生成完成')
        return hom_path


class Fom:
    # 整合一阶变异信息
    def fom_data_old(self, fompath, program_path, cov_info=[]):
        fom_record = ''
        lines = util.File().read_line(program_path)

        # 如果未给定生成变异体的行数，则默认全生成
        if not cov_info:
            cov_info = range(1, len(lines) + 1)

        # print(cov_info)
        # print(lines)
        for i in cov_info:
            line = lines[int(i) - 1]
            if self.line_pass(line):
                continue
            fom_messages = self.fom_creat(line)
            for fom_message in fom_messages:
                fom_record += str(self.fom_message_to_save([[i, line] + fom_message])) + '\n'
        if fompath == '':
            return fom_record
        with open(fompath, 'w', encoding='utf-8') as f:
            f.write(fom_record)
        return fompath

    # 整合一阶变异信息
    def fom_data(self, fompath, program_path, cov_info=[]):
        fom_record = ''
        lines = util.File().read_line(program_path)

        # 如果未给定生成变异体的行数，则默认全生成
        if not cov_info:
            cov_info = range(1, len(lines) + 1)

        # print(cov_info)
        # print(lines)
        for i in cov_info:
            line = lines[int(i) - 1]
            if self.line_pass(line):
                continue
            fom_messages = self.fomCreatNewop(line)
            for fom_message in fom_messages:
                fom_record += str(self.fom_message_to_save([[i, line] + fom_message])) + '\n'
        if fompath == '':
            return fom_record
        with open(fompath, 'w', encoding='utf-8') as f:
            f.write(fom_record)
        return fompath

    def del_nocompile_fom(self, fom_path, src_or, creat_src):
        mut_fom = ''
        chosenumlist = [0, 0, len(util.File().read_line(fom_path))]
        for j, mut in enumerate(util.File().read_line(fom_path)):
            src_mut = os.path.join(creat_src, '%s.c' % j)
            if not self.mutation_build(src_or, src_mut, eval(mut)):
                chosenumlist[1] += 1
                continue
            true_com = util.Command(src_mut)
            if not Codeflaws().mysystem(true_com.compile("%s" % j)):
                chosenumlist[1] += 1
                continue
            chosenumlist[0] += 1
            mut_fom += mut + '\n'
            # if j%20 == 0:
            #     print('选择：%s， 丢弃：%s，总共：%s' % (chosenumlist[0], chosenumlist[1], chosenumlist[2]))
        # print('选择：%s， 丢弃：%s，总共：%s' % (chosenumlist[0], chosenumlist[1], chosenumlist[2]))
        with open(fom_path, 'w') as f:
            f.write(mut_fom)
        return fom_path

    # 对一阶变异信息全随机生成二阶变异信息
    def hom_data(self, fom_path):
        print('开始生成二阶变异体')
        top = 1
        fom_record = ''
        lines = util.File().read_line(fom_path)

        toplines = {}
        inline = []
        for i, fommessage1 in enumerate(lines):
            fom1 = eval(fommessage1)[0]
            if not fom1[0] in toplines:
                toplines[fom1[0]] = 0
            if not toplines[fom1[0]] < top:
                # top限制每行取的变异体数量
                continue
            else:
                toplines[fom1[0]] += 1
                inline.append(i)

            for j, fommessage2 in enumerate(lines):
                fom2 = eval(fommessage2)[0]
                if j >= i:
                    break
                if not j in inline or fom1[0] == fom2[0]:
                    # 去除不生成组合
                    continue

                fom_record += str([fom1, fom2]) + '\n'

        with open(os.path.join(r'./report/c/Hom-original.txt'), 'w', encoding='utf-8') as f:
            f.write(fom_record)
        print('生成完成')
        return 0

    # 需要跳过不生成变异体的代码行
    def line_pass(self, line):
        passing = ['#include', 'using namespace', 'int main']
        for key in passing:
            if key in line:
                return True
        return False

    # 构造一阶变异体
    def fom_creat(self, line):
        fom_messages = []
        for key in mutation_trick[0]:
            if key in line:
                for value in mutation_trick[0][key]:
                    fom_messages.append([key, value])
        for key in mutation_trick[1]:
            if key in line:
                for value in mutation_trick[1][key]:
                    fom_messages.append([key, value])
        return fom_messages

    # 新算法构造一阶变异体
    def fomCreatNewop(self, line):
        fom_messages = []
        for operator_dict in new_mutation_trick:
            for operator_father, operator_child_list in operator_dict.items():
                if operator_father in line:
                    loc = line.find(operator_father)
                    for operator_child in operator_child_list:
                        fom_messages.append([operator_father, operator_child, loc])
                    loc = line.rfind(operator_father, loc+1)
                    while not loc == -1:
                        for operator_child in operator_child_list:
                            fom_messages.append([operator_father, operator_child, loc])
                        loc = line.rfind(operator_father, loc+1)

        return fom_messages


    # 创建存储用变异信息
    def old_fom_message_to_save(self, message):
        fom_replace = {}
        if type(message) == list:
            if '\n' == message[2][-1:]:
                before = message[2][:-1]
            else:
                before = message[2]
            after = before.replace(message[3], message[4])
            fom_replace['fom_num'] = message[0]
            fom_replace['line'] = message[1]
            fom_replace['before'] = before
            fom_replace['after'] = after
            return fom_replace
        else:
            fom_replace = eval(message)
            return [[fom_replace['fom_num'], fom_replace['line'], fom_replace['after']]]

    # 创建存储用变异信息
    def fom_message_to_save(self, messages):
        replaces = []
        for message in messages:
            single_replace = []
            if '\n' == message[1][-1:]:
                before = message[1][:-1]
            else:
                before = message[1]
            single_replace.append(message[0])
            single_replace.append(before)
            after = before[:message[4]] + message[3] + before[message[4]+len(message[2]):]
            single_replace.append(after)
            single_replace += message[2:]
            replaces.append(single_replace)
        return replaces

    # 根据变异信息进行变异
    def old_mutation_build(self, program_path, mutation_path, mutation_message):
        try:
            lines = util.File().read_line(program_path)
            lines[mutation_message['line'] - 1] = mutation_message['after']
            with open(mutation_path, 'w') as f:
                for line in lines:
                    f.write(line + '\n')
            # print("生成成功", mutation_message)
            return True
        except Exception as e:
            return False

    # 根据变异信息进行变异
    def mutation_build(self, program_path, mutation_path, mutation_message):
        try:
            # mutation_message = [[7, '     int a, p = 0, i;', '     short int a, p = 0, i;', 'int', 'short int', 5],
            #                     [7, '     int a, p = 0, i;', '     int a, p = 1, i;', '0', '1', 16]]
            lines = util.File().read_line(program_path)
            for i, om in enumerate(mutation_message):
                sameline = False
                loc = om[5]
                for j in range(0, i):
                    if mutation_message[j][0] == om[0]:
                        sameline = True
                        if mutation_message[j][5] < om[5]:
                            # 对当前index存在影响
                            loc += len(mutation_message[j][4]) - len(mutation_message[j][3])

                line_b = lines[om[0] - 1]

                if sameline:
                    lines[om[0] - 1] = lines[om[0] - 1][:loc] + om[4] + lines[om[0] - 1][loc+len(om[3]):]
                else:
                    lines[om[0] - 1] = om[2]

                line_a = lines[om[0] - 1]
            with open(mutation_path, 'w') as f:
                for line in lines:
                    f.write(line + '\n')
                # print("生成成功", mutation_message)
            return True
        except Exception as e:
            print(e)
            return False

    # 一阶变异体执行结果
    def mutation_work(self, program_path, testdata_dirpath):
        # 读取变异信息
        fom_text_original = util.File().read_line(r'./report/Fom-run.txt') if globalmessage else util.File().read_line(
            r'./report/Fom-original.txt')
        fom_text_compile = []  # 存储可通过编译的变异体信息
        fom_text_run = []  # 存储可运行的变异体信息
        fom_num = -1
        result = {}

        util.File().clear('./mbfl/fom')
        for text in fom_text_original:
            print('-----------------------------------------------------------')
            print("开始执行", text)
            fom_num += 1
            fom_message = eval(text)
            # 获取变异体执行情况
            fom_path = './mbfl/fom/Fom-%s.c' % fom_num
            util.File().remove_relative(fom_path, True)

            if not self.mutation_build(program_path, fom_path, fom_message):
                continue
            # 在fom_path位置生成相应变异体

            P2 = Programwork()
            P2.program_compile(fom_path)
            if not P2.compile:
                util.File().remove_relative(fom_path, True)
                fom_num -= 1
                continue
            fom_text_compile.append(text)
            # 对生成的变异体进行编译

            P2.program_start(fom_path, testdata_dirpath)
            if not P2.cov:
                util.File().remove_relative(fom_path, True)
                continue
            fom_text_run.append(text)
            # 对编译完成的变异体运行测试用例

            if fom_message['line'] not in result:
                result[fom_message['line']] = []
            result[fom_message['line']].append(P2.res)

            util.File().remove_relative(fom_path)
            # print("执行完成", text)

        with open(r'./report/Fom-compile.txt', 'w') as f:
            for message in fom_text_compile:
                f.write(str(message) + '\n')
        with open(r'./report/Fom-run.txt', 'w') as f:
            for message in fom_text_run:
                f.write(str(message) + '\n')
        return result


class Javafom:
    def __init__(self, src):
        self.src = os.path.join(src, 'src/main/java')
        self.creat_fom = False

    # 生成总体一阶变异信息
    def fom_data(self, cov_info):
        try:
            print('开始生成一阶变异体')
            fom_record = ''
            num = 0
            for package_name in cov_info:
                for class_path in cov_info[package_name]:
                    for line in cov_info[package_name][class_path]:
                        hits = cov_info[package_name][class_path][line]
                        if not sum(hits) > 0:
                            continue
                        filepath = os.path.join(self.src, class_path)
                        mut_operation = self.fom_info(filepath, line)
                        for before in mut_operation:
                            for after in mut_operation[before]:
                                num, message = self.fom_message(num, class_path, line, before, after)
                                fom_record += message

            with open(os.path.join(r'./report/java/Fom-original.txt'), 'w', encoding='utf-8') as f:
                f.write(fom_record)
            self.creat_fom = True
            print('创建成功')
            return True
        except:
            return False

    # 生成单个一阶变异体他信息
    def fom_info(self, path, line):
        with open(path, 'r') as f:
            text = f.readlines()[int(line) - 1]
        mut_operation = {text: []}
        for opset in java_mutation:
            for op in opset:
                slocal = text.find(op, 0)
                while not slocal == -1:
                    for other_op in opset:
                        if other_op == op:
                            continue
                        a_text = text[:slocal] + other_op + text[slocal + len(op):]
                        mut_operation[text].append(a_text)
                    slocal = text.find(op, slocal + 1)
        return mut_operation

    # 创建变异体存储格式
    def fom_message(self, fnum, fpath, fline, fbefore, fafter):
        text = str(fnum) + '----' + str(fpath) + '----' + str(fline) + \
               '----' + fbefore.strip() + '----' + fafter.strip() + '\n'
        return fnum + 1, text

    # 读取存储格式的变异体
    def read_fom_message(self, info):
        [num, path, line, before, after] = info.split('----')
        return num, path, int(line), before, after

    # 开始变异
    def startmutant(self):
        # 此时默认程序为原始程序版本

        record = util.W2py('./report/java/fom-report.txt')
        with open(r'./report/java/Fom-original.txt', 'r') as f:
            fom_oinfos = f.readlines()
        for fominfo in fom_oinfos:
            # 进行变异
            if not self.mutation_creat(fominfo, True):
                continue
            muti = util.Java('data')

            # 回退变异
            self.mutation_creat(fominfo, False)
            break

    # 根据变异信息生成变异体 可正向和反向
    def mutation_creat(self, mut_info, ty_pe):
        # ty_pe 是否正向生成变异体
        try:
            num, path, line, before, after = self.read_fom_message(mut_info)
            print(num)
            filepath = os.path.join(self.src, path)
            self.backup(filepath)
            with open(filepath, 'r') as f:
                text = f.readlines()
                if ty_pe:
                    text[line - 1] = after
                else:
                    text[line - 1] = before
                if not text[line - 1][-1] == '\n':
                    text[line - 1] += '\n'
            with open(filepath, 'w') as f:
                f.write(''.join(text))
            return True
        except:
            return False

    # 备份及恢复备份文件
    def backup(self, path):
        backpath = self.src
        # print(os.path.basename(path))
        # print(backpath)


def suspicious(mainmessage, result_fom):
    if len(mainmessage['sbfl']['testres']) == 0:
        programwork = Programwork()
        programwork.program_start(mainmessage['program_path'], mainmessage['testdata_dirpath'])
        if not programwork.cov:
            return mainmessage
        result_src = programwork.res
    else:
        result_src = mainmessage['sbfl']['testres']
    touple = get_toule(result_src, result_fom)
    sus = {}
    for key in touple['f&p']:
        # print(touple['tf'], touple['f&p'][key]['f'], touple['f&p'][key]['p'], touple['f2p'], touple['p2f'])
        sus[key] = mainmessage['mbfl']['mbflform']([touple['tf'], touple['f&p'][key]['f'], touple['f&p'][key]['p'],
                                                    touple['f2p'], touple['p2f']])
        mainmessage['mbfl']['touple'][key] = [touple['tf'], touple['f&p'][key]['f'], touple['f&p'][key]['p'],
                                              touple['f2p'], touple['p2f']]
    mainmessage['mbfl']['sus'] = sus
    return mainmessage


# 根据执行结果生成元组
def get_toule_lod(result_src, result_fom):
    # print(result_src)
    # print(result_fom)
    touple = {
        'tf': 0,
        'tp': 0,
        'f2p': 0,
        'p2f': 0,
        'f&p': {},
    }
    model = {
        'f': [],
        'p': [],
    }
    for res in result_src:
        if not res:
            touple['tf'] += 1
        else:
            touple['tp'] += 1

    for line in result_fom:
        touple['f&p'][line] = {'f': [], 'p': []}
        for m in result_fom[line]:
            f, p = 0, 0
            for i, res in enumerate(m):
                if res:
                    if not result_src[i]:
                        f += 1
                else:
                    if result_src[i]:
                        p += 1
            touple['f2p'] += f
            touple['p2f'] += p
            touple['f&p'][line]['f'].append(f)
            touple['f&p'][line]['p'].append(p)
    return touple


# 根据执行结果生成元组
def get_toule(result_or, result_src, result_fom):
    # print(result_src)
    # print(result_fom)
    touple = {
        'tf': 0,
        'tp': 0,
        'f2p': 0,
        'p2f': 0,
        'linedata': {},
    }

    for i in range(len(result_or)):
        if result_or[i] == result_src[i]:
            touple['tf'] += 1
        else:
            touple['tp'] += 1


    for line in result_fom:
        touple['linedata'][line] = {'f': [],
                                    'p': [],
                                    'kf': [],
                                    'kp': [],
                                    'nf': [],
                                    'np': []}
        for i, mut_i in enumerate(result_fom[line]):
            f, p = 0, 0
            kf, kp, nf, np = 0, 0, 0, 0
            for j, res_j in enumerate(mut_i):
                if result_or[j] == result_src[j]:
                    # 原始pass
                    if not result_or[j] == res_j:
                        # 变异体p -> f
                        f += 1
                else:
                    # 原始fail
                    if result_or[j] == res_j:
                        # 变异体 f -> p
                        p += 1

                if result_or[j] == result_src[j]:
                    # 原始pass
                    if result_src[j] == res_j:
                        # 变异体n
                        np += 1
                    else:
                        # 变异体k
                        kp += 1
                else:
                    # 原始fail
                    if result_or[j] == res_j:
                        if result_src[j] == res_j:
                            # 变异体n
                            nf += 1
                        else:
                            # 变异体k
                            kf += 1

            touple['linedata'][line]['f'].append(f)
            touple['linedata'][line]['p'].append(p)
            touple['linedata'][line]['kf'].append(kf)
            touple['linedata'][line]['kp'].append(kp)
            touple['linedata'][line]['nf'].append(nf)
            touple['linedata'][line]['np'].append(np)
        touple['f2p'] += f
        touple['p2f'] += p
    return touple


# 根据执行结果生成元组
def GetTouleList(or_list, fom_out_dic, fom_kill_dic):
    # print(result_src)
    # print(result_fom)
    touple = {
        'tf': 0,
        'tp': 0,
        'f2p': 0,
        'p2f': 0,
        'linedata': {},
    }

    for i in range(len(or_list)):
        touple['tf'] += 1 - or_list[i]
        touple['tp'] += or_list[i]

    for line in fom_out_dic:
        touple['linedata'][line] = {'f': [],
                                    'p': [],
                                    'kf': [],
                                    'kp': [],
                                    'nf': [],
                                    'np': []}
        for i in range(len(fom_out_dic[line])):
            mut_i_out = fom_out_dic[line][i]
            mut_i_kill = fom_kill_dic[line][i]
            f, p = 0, 0
            kf, kp, nf, np = 0, 0, 0, 0
            for j in range(len(or_list)):
                if or_list[j] == 1:
                    # 原始pass
                    if mut_i_out[j] == 0:
                        # 变异体p -> f
                        f += 1
                else:
                    # 原始fail
                    if mut_i_out[j] == 1:
                        # 变异体 f -> p
                        p += 1

                if or_list[j] == 1:
                    # 原始pass
                    if mut_i_kill[j] == 1:
                        # 变异体n
                        np += 1
                    else:
                        # 变异体k
                        kp += 1
                else:
                    # 原始fail
                    if mut_i_kill[j] == 1:
                        # 变异体n
                        kf += 1
                    else:
                        # 变异体k
                        nf += 1

            touple['linedata'][line]['f'].append(f)
            touple['linedata'][line]['p'].append(p)
            touple['linedata'][line]['kf'].append(kf)
            touple['linedata'][line]['kp'].append(kp)
            touple['linedata'][line]['nf'].append(nf)
            touple['linedata'][line]['np'].append(np)
            # touple['p2f'] += f
            # touple['f2p'] += p

    plist = [0 for i in range(len(or_list))]
    flist = [0 for i in range(len(or_list))]
    for line in fom_out_dic:
        for out in fom_out_dic[line]:
            for i, pf in enumerate(out):
                if pf == 0:
                    flist[i] += 1
                else:
                    plist[i] += 1

    for i in range(len(or_list)):

        if or_list[i] == 0:
            if plist[i] > 0:
                touple['f2p'] += 1
        else:
            if flist[i] > 0:
                touple['p2f'] += 1
    return touple


# 根据原始信息计算怀疑度
def getsuslvalue(touple_dic, mbfl_method, type):
    mbfl_for.type_mbfl = type
    sus_list = []
    for line in touple_dic['linedata']:
        touple_list = touple_dic['tf'], \
                      touple_dic['tp'], \
                      touple_dic['linedata'][line]['f'], \
                      touple_dic['linedata'][line]['p'], \
                      touple_dic['f2p'], \
                      touple_dic['p2f'], \
                      touple_dic['linedata'][line]['kf'], \
                      touple_dic['linedata'][line]['kp'], \
                      touple_dic['linedata'][line]['nf'], \
                      touple_dic['linedata'][line]['np'],
        sus = eval('mbfl_for.' + mbfl_method)(touple_list)
        sus_list.append([line, sus])
    return sorted(sus_list, key=lambda x: x[1], reverse=True)


