## Complex Fault Type Analysis

Here we describe the analysis of complex fault types in more detail.

### 1. Multiple Unrelated Faults on Multiple Lines

![image](https://github.com/yzzupgo/FL/blob/main/CFTA/img/Table%20I.png)

It represents that simple faults which are independent and unrelated are spread over multiple lines of a program. In such a case, failures that are observed by different test cases may be produced by different faults. However,  the results of these test cases may be affected when multiple faults are present simultaneously.

As can be seen in Table I, the faults occur on line 3 and line 10 where the statements should be ```m = z``` and ```if(x > y);``` respectively. Note that, we list the test results on each fault line when each fault exists alone and the last line represents the test result of the complex fault. The failed test cases are highlighted. When the first fault exits alone, the $t_4$ and $t_5$ are failed test cases; When the second fault exits alone, the $t_3$ is failed test case. Therefore, multiple faults are unrelated and detected by different test cases. However, the $t_4$ is passed when two faults present simultaneously.

For such type of complex fault, developers may be used to use all test cases to locate the faults if they don't distinguish the faults. However, locating a simple fault by failed test case caused by another simple fault may obtain a very poor result. Therefore, we can mutate multiple faulty statements to generate First-Order Mutants (FOMs) and combine them as a Higher-Order Mutant (HOM). The HOM can simulate such type of complex fault if it can pass all test cases. In summary, HOM is a potential way to fix complex fault. 

### 2. Multiple Related Faults on Multiple Lines

![image](https://github.com/yzzupgo/FL/blob/main/CFTA/img/Table%20II.png)

It represents that simple faults which are related are spread over non-contiguous statements of a program. Note that, the statements may distribute in different functions and different classes of C++ and Java.

As can be seen in Table II, the faults occur on line 4 and line 6 where the statements should be ```if(y < z)``` and ```m = y;``` respectively. Specially, the multiple faults occur under an ```if``` branch, which determines the faults are related. Furthermore, when the first fault exits alone, the $t_2$ and $t_5$ are failed test cases; When the second fault exits alone, the $t_2$ is a failed test case. Therefore, all test cases can pass only when multiple faults are fixed.

For such type of complex fault, developers can completely fix it only when they locate and fix total simple faults. Although multiple simple faults are related, we can fix them in a similar way to the first method, i.e., generating HOM by mutating faulty statements to simulate such type of complex fault.

### 3. Multiple Related Faults on A Single Line

![image](https://github.com/yzzupgo/FL/blob/main/CFTA/img/Table%20III.png)

It represents that there are multiple faults on a single line with a certain fault relevance.

As can be seen in Table III, the faults only occur on line 4 where the statement should be ```if(y < z)```. Multiple faults are ```<=``` and ```z - 1``` respectively. When the first fault exits alone, the $t_1$, $t_2$, $t_3$, $t_5$ and $t_6$ are failed test cases; When the second fault exits alone, the $t_2$ and $t_5$ are failed test cases. However, when the multiple faults exit together, the $t_1$, $t_3$ and $t_6$ are failed test cases while $t_2$ and $t_5$ change into passed test cases. Therefore, the multiple faults are related and test cases can't distinguish them.

For such type of complex fault, failed test cases can't distinguish multiple simple faults as the faults occur on a single line. Furthermore, the faults seem to be caused by the same cause from their performance. To fix such type of complex fault, we can apply HOM on a single line, i.e., mutating the faulty statement multiple times.


### Conclusion

Based on the above studies, we can draw that the complex fault consisting of multiple simple faults and analyzing the relationship between simple faults to generate an optimal HOM is an effective way to simulate complex fault. According to the relevance and distribution characteristic of simple faults, we propose three HOM generation methods (i.e., SFClu, SFDis, and SFDen).
