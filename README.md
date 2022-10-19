# lunwen

## Complex Fault Type Analysis

Here we describe the analysis of complex fault types in more detail.

### 1. Multiple Unrelated Faults on Multiple Lines

![image](https://github.com/yzzupgo/lunwen/blob/main/img/Table%20I.png)

For such type of complex fault, developers may be used to use all test cases to locate the faults if they don't distinguish the faults. However, locating a simple fault by failed test case caused by another simple fault may obtain a very poor result. Therefore, we can mutate multiple faulty statements to generate First-Order Mutants (FOMs) and combine them as a Higher-Order Mutant (HOM). The HOM can simulate such type of complex fault if it can pass all test cases. In summary, HOM is a potential way to fix the complex fault. 

### 2. Multiple Related Faults on Multiple Lines

![image](https://github.com/yzzupgo/lunwen/blob/main/img/Table%20II.png)

For such type of complex fault, developers can completely fix it only when they locate and fix total simple faults. Although multiple simple faults are related, we can fix them in a similar way to the first method, i.e., generating HOM by mutating faulty statements to simulate such type of complex fault.

### 3. Multiple Related Faults on A Single Line

![image](https://github.com/yzzupgo/lunwen/blob/main/img/Table%20III.png)

For such type of complex fault, failed test cases can't distinguish multiple simple faults as the faults occur on a single line. Furthermore, the faults seem to be caused by the same reason from their performance. To fix such type of complex fault, we can apply HOM on a single line, i.e., mutating the faulty statement multiple times.


### Conclusion

Based on the above studies, we can draw that the complex fault consisting of multiple simple faults and analyzing the relationship between simple faults to generate an optimal HOM is an effective way to simulate complex fault. According to the relevance and distribution characteristic of simple faults, we propose three HOM generation methods (i.e., SFClu, SFDis, and SFDen).
