本程序旨在通过代码进行行星轮系的配齿设计，采用matlab。（English version below）
代码里展现了进行运算的具体函数，采用轮系名称进行的命名。关于行星轮系简介可见：
其中NGW作为最常见的轮系类型，直接采用暴力的历遍算法。
对于NW, NN ,WW, 这三类具有双联行星轮的“两排”轮系，由于涉及传动比分配问题，历遍算法显得非常鸡肋。
设计手册给出了NGW的传动比分配公式，提出“应以各级之间获得等强度和最小外廓尺寸为原则。”
此外也有相关文献针对轮系特别是多档位变速器传动比分配进行了研究。
不过，在设计程序时考虑到应该得出尽量多的齿数组合供用户选择，在设计过程中扩大了一定传动比分配的范围。
参考文献
[1] 张展, 张国瑞, 张焕武, 等. 实用齿轮设计手册[M]. 北京: 机械工业出版社, 2010.
[2] 饶振纲. 行星齿轮传动设计[M]. 北京: 化学工业出版社, 2014.
This program is designed to design the gearing of planetary gear trains through code, using Matlab.
The code shows the specific functions for the calculations.Naming using the name of the gear train. For an introduction to planetary gear trains, see:
NGW is the most common type of gear train, and a brute force traversal algorithm is used directly.
For the three types of "two-row" gear trains with double-linked planetary gears, NW, NN, and WW, the traversal algorithm is very useless because of the transmission ratio distribution problem.
The design manual gives the transmission ratio distribution formula of NGW, and proposes that "the principle should be to obtain equal strength and minimum external dimensions between each level."
In addition, there are also relevant literatures that have studied the transmission ratio distribution of gear trains, especially multi-speed transmissions.
However, when designing the program, considering that as many tooth number combinations as possible should be obtained for users to choose from, 
the range of certain transmission ratio distribution was expanded during the design process.
References
[1] Zhang Zhan, Zhang Guorui, Zhang Huanwu, et al. Practical Gear Design Manual [M]. Beijing: Machinery Industry Press, 2010.
[2] Rao Zhengang. Planetary Gear Transmission Design [M]. Beijing: Chemical Industry Press, 2014.
