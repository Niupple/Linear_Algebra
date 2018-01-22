# 数学库的设计：Linear\_Algebra 
## 概况 
### 库介绍

本数学库提供了若干符合C++使用习惯的代数结构及其常用计算。

本数学库是Yuanjie "Duane" Ding @Niupple的一个小项目，bug可能并没有de完，函数设计也可能并不尽人意，但是如果真的有人要用这个库的话，此帮助文档可以帮您快速地了解这个库的作用和使用方式。此外，帮助文档包含了库中全部函数的具体实现细节，此用于个人备忘的同时，也可以帮助您全面地认识这个库的效率、局限以及可能的漏洞。

email: niupple@gmail.com

### 整体设计  

整个数学库都是由不同的、平行的代数结构组成的，它们包括抽象的：
1.	群
1.	环
1.	域

这些代数结构均为模板类，按照其惯有用法定义了乘法、加法、除法、逆元等概念，并且按照数学上对这些结构的概念大小，规定了它们的基类和一系列派生类。
同时，我们还具体实现了一些经常研究的、有重要意义的代数结构，例如：
1.	整数环（无限精度的）`<integer.h>`
1.	有理数域（基于整数环的）`<rational.h>`
1.	一般线性群（基于某个数域的）`<linear_algebra.h>`
1.	多项式环（在某数域上的）`<polinomial.h>`

并且针对这些代数结构的独有性质，我们实现了诸如整数环上的同余运算、在模互质数意义下的逆元运算、有理数域上的倒数运算、模有理数、方阵求行列式、可逆矩阵求逆等等一系列运算。

在代码风格方面，我们沿用了我一贯的风格，并且在函数、类命名上采用了C++风格的下划线分隔命名法，翻译一概采用自维基百科的英文页面，并在本文档中给出了翻译对照表。为了编辑的方便，代码中的注释以及说明文本一概使用*英文书写*。

在构造方面，我们期望这些模板类使用起来就像内置类型一样便利，于是我们在设计中尽量充分地考虑了使用时可能需要用到的构造函数、运算符重载，使得它们的使用如同内置类型一样灵活而定义良好。

在输入输出方面，我们期望诸如**整数**、**有理数**、**多项式**这样的，与内置数据类型相似的类能够适配`cin`、`cout`等流控制的输入输出方式；而**矩阵**等与标准库容器相似的类适配迭代器，而不整体适配输入输出流。

在函数封装方面，我们封装了很多非常好用的函数，基本可以使用最少量的参数就进行完整的运算。为了进一步增加模板类使用的灵活性，我们为每一个可能修改类的函数和成员函数提供了针对左值的修改版本和针对右值的拷贝版本，并且力求用户完全知晓并自然地注意到他们正在进行的操作是一个修改参数的，还是一个返回拷贝的。

在成员可访问性方面，我们坚持了使用重载的下标运算符`[idx]`和特定的成员函数的间接访问，类内保存的核心数据及类维护的信息是不能够在外部访问的。

## 作用

### 整型`<intger.h>`

### 有理数`<rational.h>`

`<rational.h>`提供了C++并不支持的，精确的有理数类，它的所有相关的类、函数的实现都在命名空间`rat`中。使用`rat::rational`类，可以进行基本的构造、拷贝、输入、输出、数学运算，也可以进行诸多的类型转换。

#### 构造函数

`rat::rational::rational()`
是**默认的**默认构造函数，它会生成一个值为0的有理数。

```cpp
cout << rational::rational() << endl;
```

`rat::rational::rational(const rat::rational&)`
是**默认的**拷贝构造函数，它会将形参完全地拷贝。

```cpp
rat::rational ra;
cin >> ra;
rat::rational rb(ra);   //rb将会与ra拥有完全相同的属性
```

`rat::rational::rational(const hint &a, const hint &b = 1)`
是根据两个`hint`类型的整型为分子（第一个形参）和分母，构造有理数的构造函数，其中由于第二个形参是带有默认值1的，它同时可以将某一个整型转换为分母为1的有理数。
函数并不限制两个参数中是否有负数，如果存在的话，函数仍然将返回值定义为 $$x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}$$ 

```cpp
rat::rational   ra(3),                  //正确，ra的值将等于3
                rb(3, 2), rc(-3, -2),   //正确，rb，rc的值都是3/2
                rd(-3, 2), re(3, -2),   //正确，rd，re的值都是-3/2
                rf(0, 3), rg(0, -1),    //正确，rf，rg的值都是0
                rh(0, 0), ri(1, 0);     //错误，第二个参数绝对不能为零！
```

#### 拷贝

`rational`类的拷贝行为是**默认的**，它会把右值的信息完全地复制到左值中。

#### 输入、输出


