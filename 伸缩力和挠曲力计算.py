import math
import sys
import matplotlib.pyplot as plt
sys.setrecursionlimit(900)  # 这里设置大一些
'''
#  程序名:轨道强度检算
#  作  者:李建英
#         中南大学土木工程学院
#          2020年5月10日02:00:39
#
  说   明:铁路轨道伸缩力和挠曲力计算
#  功能说明：
#  修改备注:
'''
class constant0():
    # 用于数据，坚守函数间传递的变量数量
    pass

CC = constant0()

# l0的编辑计算过程
CC.num_arc = 7
CC.r_Rail = 43.3 # 桥两端阻力
CC.r0 =70.0  # 桥上线路阻力
CC.r1 = 70
CC.r2 = 66.7
CC.r3 = 43.3
CC.r4 = 43.3
CC.l1 = 750   #机车长度

CC.brige_arclenth = [3200 for num in range(CC.num_arc)]
#CC.brige_arclenth = [3200 for num in range(CC.num_arc-1)]  # 桥跨长度
#CC.brige_arclenth.append(2400)
CC.left_type = [0 for num in range(CC.num_arc)]  # 需要左端支座类型外部数据，来进行判断（固定为0，可动为1）
CC.brige_flex = { 0:0 , 200:0.011 , 400:0.042 , 600:0.091 , 800:0.154 , 1000:0.228 , 1200:0.311 , 1400:0.400 , 1600:0.492 , 1800:0.548 ,2000:0.673 , 2200:0.756 , 2400:0.830 , 2600:0.893 , 2800:0.942 ,3000:0.973 ,3200:0.984}
#说明：距活动端距离（cm）：位移（cm）

CC.alpha =1.0e-5
CC.delta_T = 15.0
CC.E_rail = 2.1e5
CC.F_rail = 77.45

def initJudgeCalculation(pCC):
    global CC
    CC = pCC    #存储全局需要的变量



def contractility(l0):  # 计算伸缩力的部分
    global l_num
    global P_num
    global ya_num
    global yb_num    #用于画图，每次计算都将记录
    temp = 0    #寻找最小值的过程中用来临时存储符合要求的l0的变量
    kesai =0.01     #工程要求闭合精度
    rule = 0.1        #计算步长，每次都加一个rule进行暴力搜索
    temp_delta = 0   #寻找最小值的过程中用来临时存储中间结果闭合精度的变量，非常大的数值表示一定不合适
    flag = 0        #计算出结果的标志，如果flag是0，表示搜索过程中并没有合适值，避免未定义的参数调用错误（很多参数，都是要符合要求后才会赋值）

    for i in range(100000):   #在l0到range范围内暴力搜索
        l_num = [0.0]
        P_num = [0.0]
        ya_num = [0.0]
        yb_num = [0.0]
        delta_y = jumponarc(l0)     #初始化储存画图用的坐标点列表

        #if math.fabs(delta_y)<kesai and math.fabs(delta_y)<temp_delta:
        if math.fabs(delta_y) < kesai and l0 > temp_delta:
            temp_delta=l0
            temp = l0
            l_last = [0.0]
            P_last = [0.0]
            ya_last = [0.0]
            yb_last = [0.0]     #初始化最终使用的画图的坐标列表
            flag+=1         #表示已经进入符合精度的要求
            for j in range(len(l_num)-1):
                l_last.append (l_num[j+1])
                P_last.append (P_num[j+1])
                ya_last.append (ya_num[j+1])
                yb_last.append (yb_num[j+1])     #如果满足要求的话，就把绘图要求的点赋值到最后使用的画图列表

        #print(l0, delta[i])
        l0 = l0+rule
    if flag!=0:
        l0 = temp
        print(l0)
        print('最大值'+str(max(P_last))+'\n'+'最小值'+str(min(P_last))+'\n')
        l_last = num_recover(l0, l_last)
        paint(l_last, P_last,ya_last,yb_last)

    else:
        print("None")


def jumponarc(l0):  # l0为路基变形长度，CC.num_arc为桥跨数
    global y,a_front,b_front,Sk1,P1,lk,Pk,flag_ym  #用于传递计算时需要的数据

    global l_num,P_num,ya_num,yb_num    #用于画图，每次计算都将记录

    P1 = CC.r0 * l0  # 表面上是P1的大小，其实接下来是作为桥跨左端点的力值不断地传递
    Sk1 = (CC.r0*l0**2)/2.0  # 桥左端的路基阻力图面积

    Pk = [0 for num in range(CC.num_arc)]
    lk = [0 for num in range(CC.num_arc)]
    b_front = Sk1/( CC.E_rail * 100 * CC.F_rail)  # b是钢轨位移，为了节省输入，用b代替，b为全局变量，每次修改都会累计

    for index in range(len(CC.left_type)):  # 在整个桥上进行循环计算
        if CC.left_type[index] == 0:  # 简支梁左端务必在每次循环前，将σa初始化
            a_front = 0.0  # a是粱位移，每次循环都会重置
            CC.r = CC.r_Rail

            l_num.append(0)
            P_num.append(P1)
            yb_num.append(b_front)
            ya_num.append(a_front)

            judge0(0,index)
        elif CC.left_type[index] == 1:
            a_front = CC.alpha * CC.delta_T * CC.brige_arclenth[index]  # a是粱位移，每次循环都会重置
            if a_front > b_front: #粱位移大于钢轨位移时，阻力方向向右，符号变为负
                CC.r = 0 - CC.r_Rail
            else:
                CC.r = CC.r_Rail

            l_num.append(0)
            P_num.append(P1)
            yb_num.append(b_front)
            ya_num.append(a_front)

            judge1(0,index)
    l_num.append(CC.brige_arclenth[index]+math.fabs(P1/CC.r_Rail))
    P_num.append(0)
    yb_num.append(0)
    ya_num.append(0)

    y = b_front + ((P1 *math.fabs(P1)) / (2 * CC.E_rail * 100 * CC.F_rail *CC.r_Rail))
    if flag_ym == 1:
        return y
    else :
        y = 10000
        return y

def judge0(lx, index):  # lx表示当前计算的位置
    global a_front, b_front,Sk1,P1,Pk,y,flag_ym
    global l_num, P_num, ya_num, yb_num  # 用于画图，每次计算都将记录
    flag_ym = 0
    while lx < (CC.brige_arclenth[index] - 200.0):
        a_rear = a_front
        b_rear = b_front
        lx += 200.0
        a_front = CC.alpha * CC.delta_T * lx
        b_front = (
            Sk1 + (2 * P1 + CC.r * (lx - lk[index]))*(lx - lk[index])/2) / (CC.E_rail * 100 * CC.F_rail)
        if (a_front - b_front) * (a_rear - b_rear) < 0:   # 判断变号
            '此时联立方程组α*δT*lx=（Sk1+(2P1+rlx）*lx/2)/EF'
            '整理方程，得到：(r/2)a²+(P1-EFα*δT)a+Sk1=0'
            delta_lx = (P1-CC.E_rail * 100 * CC.F_rail *
                        CC.alpha*CC.delta_T)**2.0-4*(CC.r/2)*Sk1
            flag_ym+=1
            if delta_lx < 0:
                y=10000
                break
            elif delta_lx > 0:
                lx1 = (-(P1-CC.E_rail * 100 * CC.F_rail*CC.alpha *
                         CC.delta_T) + math.sqrt(delta_lx)) / (2 * CC.r/2)
                lx2 = (-(P1-CC.E_rail * 100 * CC.F_rail*CC.alpha *
                         CC.delta_T) - math.sqrt(delta_lx)) / (2 * CC.r/2)
                if lx1 > 0 and lx1 < CC.brige_arclenth[index]:
                    if lx2 > 0:
                        lk[index] = min(lx1, lx2)
                    else:
                        lk[index] = lx1
                elif lx2 > 0 and lx2 < CC.brige_arclenth[index]:
                    if lx1 > 0:
                        lk[index] = min(lx1, lx2)
                    else:
                        lk[index] = lx2
                else:
                    y=10000
                    break
            elif delta_lx == 0:
                lk[index] = (-(P1 - CC.E_rail * 100 * CC.F_rail *
                               CC.alpha * CC.delta_T)) / (2 * CC.r / 2)
            Sk1 = Sk1+(2*P1+CC.r*lk[index])*lk[index]/2
            P1 = P1 + CC.r * lk[index]
            Pk[index] = P1
            CC.r = 0-CC.r
            a_front = CC.alpha * CC.delta_T * lk[index]
            b_front =a_front
            lx=lk[index]
        #结束位移相等点的处理
        l_num.append(lx)
        P_num.append(P1+CC.r*(lx-lk[index]))
        ya_num.append(a_front)
        yb_num.append(b_front)
    else:                                                                                                                                                                           # coding by jianying li

        Sk1 = Sk1 + (2 * P1 + CC.r * (CC.brige_arclenth[index]  - lk[index])) * \
            (CC.brige_arclenth[index] - lk[index])/2
        b_front = Sk1 / (CC.E_rail * 100 * CC.F_rail)
        a_front = CC.alpha * CC.delta_T*CC.brige_arclenth[index]
        P1 = P1 + CC.r * (CC.brige_arclenth[index] - lk[index])
        l_num.append(CC.brige_arclenth[index])
        P_num.append(P1)
        ya_num.append(a_front)
        yb_num.append(b_front)

def judge1(lx,index):
    global a_front,b_front,Sk1,P1,lk,Pk,y
    global l_num, P_num, ya_num, yb_num  # 用于画图，每次计算都将记录
    while lx < (CC.brige_arclenth[index] - 200.0):
        a_rear = a_front
        b_rear = b_front
        lx += 200.0
        a_front = CC.alpha * CC.delta_T * (CC.brige_arclenth[index] - lx)
        b_front = (
            Sk1 + (2 * P1 + CC.r * (lx - lk[index])*(lx - lk[index])/2) / CC.E_rail * 100 * CC.F_rail)
        if (a_front - b_front) * (a_rear - b_rear) < 0:
            # '此时联立方程组α*δT*lx=（Sk1+(2P1+rlx）*lx/2)/EF'
            # '整理方程，得到：(r/2)lx²+(P1+EFα*δT)lx+Sk1-EFα*δTlen=0'
            delta_lx = (P1 + CC.E_rail * 100 * CC.F_rail * CC.alpha * CC.delta_T) ** 2.0 - \
                4 * (CC.r / 2) * (Sk1 - CC.E_rail *
                                  CC.F_rail * CC.alpha * CC.brige_arclenth[index])
            if delta_lx < 0:
                y=10000
                break
            elif delta_lx > 0:
                lx1 = (-(P1+CC.E_rail * 100 * CC.F_rail*CC.alpha *
                         CC.delta_T) + math.sqrt(delta_lx)) / (2 * CC.r/2)
                lx2 = (-(P1+CC.E_rail * 100 * CC.F_rail*CC.alpha *
                         CC.delta_T) - math.sqrt(delta_lx)) / (2 * CC.r/2)
                if lx1 > 0 and lx1 < CC.brige_arclenth:
                    if lx2 > 0:
                        lk[index] = min(lx1, lx2)
                    else:
                        lk[index] = lx1
                elif lx2 > 0 and lx2 < CC.brige_arclenth:
                    if lx1 > 0:
                        lk[index] = min(lx1, lx2)
                    else:
                        lk[index] = lx2
                else:
                    y=10000
                    break
            elif delta_lx == 0:
                lk[index] = (-(P1+CC.E_rail * 100 * CC.F_rail *
                               CC.alpha * CC.delta_T)) / (2 * CC.r / 2)

            Sk1 = Sk1+(2*P1+CC.r*lk[index])*lk[index]/2
            CC.r = 0-CC.r
            P1 = P1 + CC.r * lk[index]
            Pk[index] = P1
            a_front = CC.alpha * CC.delta_T * \
                (CC.brige_arclenth[index] - lk[index])
            b_front = Sk1 / (CC.E_rail * 100 * CC.F_rail)
            lx = lk[index]
        l_num.append(lx)
        P_num.append(P1+CC.r*(lx-lk[index]))
        ya_num.append(a_front)
        yb_num.append(b_front)
    else:

        Sk1 = Sk1 + (2 * P1 + CC.r * (CC.brige_arclenth[index] - lk[index])) * \
            (CC.brige_arclenth[index] - lk[index])/2
        b_front = Sk1 / (CC.E_rail * 100 * CC.F_rail)
        a_front = 0
        P1 = P1 + CC.r * (CC.brige_arclenth[index] - lk[index])
        l_num.append(CC.brige_arclenth[index])
        P_num.append(P1)
        ya_num.append(a_front)
        yb_num.append(b_front)

def num_recover(l0,num):
    temp = 0
    for index in range(1,len(num)):   #第一个数是0，所以从第二个数，即桥跨起点开始算0点
        if num[index] == 0:
            temp = num[index-1]
        num[index] +=temp

    num[0] = 0-l0
    return num

#挠曲力计算说明
def flexing_force(l0):    #计算挠曲力的部分
    global x_num, Px_num, yr_num,len_brige
    kesai = 0.001
    rule = 1
    temp_delta = 1000
    flag = 0
    for i in range(4000):

        delta_y = flexing_judge(l0)
        if math.fabs(delta_y) < kesai and math.fabs(delta_y) < temp_delta:
            temp_delta = math.fabs(delta_y)   #保存最小值
            temp = l0
            x_last = [x_num[0]]
            Px_last = [0.0]
            yr_last = [0.0]
            flag+=1
            for j in range(len(x_num) - 1):
                x_last.append(x_num[j + 1])
                Px_last.append(Px_num[j + 1])
                yr_last.append(yr_num[j + 1])

            # print(l0, delta[i])
        l0 = l0 + rule
    if flag != 0:
        l0 = temp
        print(str(l0)+'\n')
        paint1(x_last,Px_last,yr_last,len_brige)

    else:
        print("None")


def flexing_judge(l0):
    global P1,P2,y1,y2,len_brige
    global x_num, Px_num, yr_num
    #len_brige = max(CC.brige_arclenth)
    len_brige = 3200
    P1 = CC.r1*l0
    y1 = P1*l0/(2*CC.E_rail * 100 * CC.F_rail)
    P2 = P1-CC.r2*CC.l1
    y2 = (P1+P2)*CC.l1/(2*CC.E_rail * 100 * CC.F_rail)  #钢轨位移 单位cm
    l1_left = (CC.l1/100)//2 *200   #向左方向取整距离单位长、cm
    yl1_left = CC.brige_flex[len_brige - l1_left]   #梁位移 单位cm
    yi = CC.brige_flex[len_brige-200 - l1_left]
    yl1 = yl1_left + ((yl1_left - yi)*(l1_left-CC.l1))/200  #粱位移 单位cm

    if y1+y2 < yl1:    #说明：计算书上用i点（距左端8米处）梁位移与钢轨位移比较，这里使用l1处粱位移（线性内插法）与钢轨位移比较
        for li in range(len_brige-200-int(l1_left),0,-200):
            if K_R(li)*K_R(li-200) <0:  #表示存在K点
                x = len_brige-CC.l1
                y = y1+y2
                z = (CC.brige_flex[li]-CC.brige_flex[(li-200)])/200
                '根据位移协调方程整理得：r3*lk**2 + 2*(EFz+P2-r3*x)*lk + 2EF(δi-zli-y)+r3*x**2-2P2*x = 0'
                a = CC.r3
                b = 2*(CC.E_rail * 100 * CC.F_rail*z + P2 - CC.r3*x)
                c = 2*CC.E_rail * 100 * CC.F_rail*(CC.brige_flex[li] - z*li -y)+CC.r3*x*x-2*P2*x
                delta_lk = b**2 - 4*a*c
                if delta_lk < 0:
                    return 10000
                    break
                elif delta_lk > 0:
                    lk1 = (-b + math.sqrt(delta_lk)) / (2 *a)
                    lk2 = (-b - math.sqrt(delta_lk)) / (2 *a)
                    if lk1 > 0 and lk1 < len_brige-200:
                        if lk2 > 0:
                            lk = min(lk1, lk2)
                            return tools(lk,l0)
                            break
                        else:
                            lk = lk1
                            return tools(lk,l0)
                            break
                    elif lk2 > 0 and lk2 < len_brige-200:
                        if lk1 > 0:
                            lk= min(lk1, lk2)
                            return tools(lk,l0)
                            break
                        else:
                            lk = lk2
                            return tools(lk,l0)
                            break
                    else:
                        return 10000
                        break
                elif delta_lk == 0:
                    lk = -b/(2*a)
                    return tools(lk,l0)
                    break
        else:
            return 10000

def tools(lk,l0):
    global P1, P2, y1, y2, len_brige
    global x_num,Px_num,yr_num
    Pk = P2 - CC.r3*(len_brige-CC.l1-lk)
    yk = y1+y2 + ((P2+Pk)*(len_brige-CC.l1-lk)/(2*CC.E_rail * 100 * CC.F_rail))
    P3 = Pk+CC.r3 *lk
    y3 = (Pk+P3)*lk/(2*CC.E_rail * 100 * CC.F_rail)
    y4 = P3*math.fabs(P3)/(2*CC.E_rail * 100 * CC.F_rail*CC.r4)

    #计算该点生成的挠曲力和钢轨纵向位移
    x_num=[0.0]
    x_num[0]=0-l0
    Px_num=[0.0]
    yr_num=[0.0]
    i=0
    while i < CC.l1 :
        x_num.append(i)
        Pi = P1 - CC.r2 *i
        Px_num.append(Pi)
        yi = y1+(P1 +Pi)*i/(2*CC.E_rail * 100 * CC.F_rail)
        yr_num.append(yi)
        i+=200
    else :
        x_num.append(CC.l1)
        Px_num.append(P2)
        yr_num.append(y1 + y2)


    while i < len_brige - lk :
        x_num.append(i)
        Pi = P2 - CC.r3 *(i-CC.l1)
        Px_num.append(Pi)
        yi = y1 +y2+ (P2 + Pi) * (i-CC.l1) / (2 * CC.E_rail * 100 * CC.F_rail)
        yr_num.append(yi)
        i += 200
    else:
        x_num.append(len_brige - lk)
        Px_num.append(Pk)
        yr_num.append(yk)

    while i < len_brige :
        x_num.append(i)
        Pi = Pk +CC.r3 *(lk-len_brige+i)
        Px_num.append(Pi)
        yi = yk + (Pk + Pi) * (lk-len_brige+i) / (2 * CC.E_rail * 100 * CC.F_rail)
        yr_num.append(yi)
        i += 200
    else:
        x = len_brige-P3/CC.r4
        x_num.append(x)
        Px_num.append(0)
        yr_num.append(yk +y3 +y4)

    return yk +y3 +y4


def K_R(li):
    global P2,y1,y2,len_brige
    yi = y1 + y2 +((P2*2-CC.r3*(len_brige-CC.l1-li))*(len_brige-CC.l1-li)/(2*CC.E_rail * 100 * CC.F_rail))

    delta = CC.brige_flex[li]-yi
    return delta



def paint(x0,p0, y0,y1):
    '''画钢轨和梁的位移图，钢轨力图'''
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.pylab as pl

    #转换单位
    p0=[a/1000.0 for a in p0] #kN
    y0=[a*1.0 for a in y0] #mm
    y1=[a*1.0 for a in y1] #mm

    x0=[a*0.01 for a in x0] #m
    # simhei.ttf 是电脑控制面板里字体的一种，这里是黑体
    # simkai.ttf 是电脑控制面板里字体的一种，这里是楷体
    chinese = matplotlib.font_manager.FontProperties(fname='C:\Windows\Fonts\simhei.ttf')

    # 参见https://blog.csdn.net/sunshihua12829/article/details/52786144/
    plt.figure(21)  # 整个绘图区域划分成2行 1列
    plt.subplot(211)  # 本图位于2行 1列的第1行
    plt.title('伸缩力钢轨和桥梁位移图(mm)', fontproperties=chinese)
    #画钢轨位移
    pl.plot(x0, y1, 'r', label='钢轨位移（mm）')

    # 画梁位移
    y=[0.0]
    temp = 0.0
    x=[temp]
    for i in range(1,len(y0)) :
        if y0[i]==0:
            pl.plot(x, y, 'b')
            y = [0.0]
            x = [x0[i]]
        else:
            x.append(x0[i])
            y.append(y0[i])
    pl.plot(x, y, 'b',label='梁位移（mm）')
    #坐标轴调整优化
    top1 = 1.1*max(max(y0), max(y1))
    bottom1 =1.1* min(min(y0), min(y1),-0.05)
    pl.ylim(top1,bottom1)
    ax = plt.gca()
    ax.spines['bottom'].set_position(('data', 0))
    ax.spines['left'].set_position(('data', 0))
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    plt.legend(prop=chinese)
    plt.grid()
    #钢轨力
    plt.subplot(212)  # 2行 1列的第2行
    plt.title('伸缩力钢轨力图(kN)', fontproperties=chinese)
    pl.plot(x0, p0, 'r', label='钢轨力（kN）')

    for j in range(len(x0)):
        if x0[j] == 0 or x0[j] %32==0:
            print('支座处伸缩力'+str(p0[j]) + '\n')

    plt.subplots_adjust(left=0.05,bottom=0.1,right=0.95,top=0.90,wspace=None,hspace=0.19)
    plt.legend(prop=chinese)
    ax = plt.gca()
    ax.spines['bottom'].set_position(('data', 0))
    ax.spines['left'].set_position(('data', 0))
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    plt.grid()
    plt.show()

def paint1(x0,p0,y1,la):
    #挠曲力绘图部分
    '''画钢轨和梁的位移图，钢轨力图'''
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.pylab as pl

    # 转换单位
    p0 = [a / 1000.0 for a in p0]  # kN
    y1 = [a * 1.0 for a in y1]  # mm
    x0 = [a * 0.01 for a in x0]  # m
    la = la/100
    # simhei.ttf 是电脑控制面板里字体的一种，这里是黑体
    # simkai.ttf 是电脑控制面板里字体的一种，这里是楷体
    chinese = matplotlib.font_manager.FontProperties(fname='C:\Windows\Fonts\simhei.ttf')

    # 参见https://blog.csdn.net/sunshihua12829/article/details/52786144/

    plt.figure(21)  # 整个绘图区域划分成2行 1列
    plt.subplot(211)  # 本图位于2行 1列的第1行
    plt.title('挠曲力钢轨和桥梁位移图(mm)', fontproperties=chinese)

    # 画钢轨位移
    pl.plot(x0, y1, 'r', label='钢轨位移（mm）')
    # 画梁位移
    y = []
    x = []
    for i in range(len(CC.brige_flex)):
            x.append(la-2*i)
            y.append(CC.brige_flex[200*i])
    pl.plot(x, y, 'b', label='梁位移（mm）')

    top1 = 1.1*max(max(y),max(y1))
    bottom1 = 1.1*min(min(y),min(y1),-0.05)
    pl.ylim(top1, bottom1) #控制坐标轴方向，只需要把上下限界反向输入即可

    ax = plt.gca()
    ax.spines['bottom'].set_position(('data', 0))
    ax.spines['left'].set_position(('data', 0))
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')

    plt.legend(prop=chinese)
    plt.grid()
    # 钢轨力
    plt.subplot(212)  # 2行 1列的第2行
    plt.title('挠曲力钢轨力图(kN)', fontproperties=chinese)
    pl.plot(x0, p0, 'r', label='钢轨力（kN）')

    #for j in range(len(x0)):
        #if x0[j] == 0:
            #print('支座处挠曲力'+str(p0[j]) + '\n')
            #print('支座处挠曲力+6m' + str(p0[j+3]) + '\n')

    ax = plt.gca()
    ax.spines['bottom'].set_position(('data', 0))
    ax.spines['left'].set_position(('data', 0))
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')

    plt.subplots_adjust(left=0.05,bottom=0.1,right=0.95,top=0.90,wspace=None,hspace=0.19)
    plt.legend(prop=chinese)
    plt.grid()
    plt.show()


if __name__ == '__main__':
    #contractility(0)
    flexing_force(0)