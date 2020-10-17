import os,sys,ogr,gdal,hashlib
from math import *
from time import *
from Crypto.Cipher import AES
from binascii import b2a_hex, a2b_hex
import numpy as np

os.chdir(r'D:\2020\20200824 加密算法英文稿子\加密算法英文 稿子、数据与程序\test data')

def gene_hashkey(user_key):
    # 密钥K的生成
    sha = hashlib.sha512()
    sha.update(user_key.encode('utf-8'))
    hash_key = sha.hexdigest()
    value_list = [int(i, 16) for i in list(hash_key)]
    return value_list

def Get_XY_fromshp(ori_shp):

    X, Y = [], []
    a = []
    ds = ogr.Open(ori_shp,0)
    in_layer = ds.GetLayer(0)
    feature_num = in_layer.GetFeatureCount()

    for in_feature in in_layer:

        geom = in_feature.geometry()
        wkt = geom.ExportToWkt()

        pointstring = wkt[wkt.rfind('(') + 1:wkt.find(')')].split(',')

        x, y = [],[]
        for point_value in pointstring:
            if point_value == '':
                continue
            else:
                xy = point_value.split(' ')

                x += [float(xy[0])]
                y += [float(xy[1])]

        X.append(x), Y.append(y)

        for i in range(1, in_feature.GetFieldCount()):
            value = in_feature.GetField(i)
            a += [value]
    # b = np.array(a).reshape(1241,9)#b是属性值矩阵化
    # fclass = b[...,1]
    # name = b[...,2]
    # print(name)
        # print(wkt.find('(98.4127828'))
    # print(wkt[wkt.find('((') + 1:wkt.rfind(')')].split(','))
    # print(len(X),feature_num)
    # print(X[0],X[-1])
    # print(Y[0],Y[-1])
    ds.Release()
    return X,Y,feature_num

def encryption_XY(ori_shp,user_key):

    begin_time = time()

    hash_key = gene_hashkey(user_key)
    X,Y,feature_num = Get_XY_fromshp(ori_shp)
    En_X, En_Y = [],[]
    # mode = AES.MODE_ECB
    # cryptos = AES.new(key, mode)
    for feat_count in range(feature_num):

        A, B = [],[]
        en_x, en_y = [],[]
        ran_x, ran_y = [],[]
        # DFT,IDFT = 0,0
        # idft_Al, idft_Alx, idft_Aly = [], [], []
        for point_num in range(len(X[feat_count])):
            a = ((feat_count+1) * (point_num+1)) / feature_num
            b = 1 / sqrt(2 * pi)
            if point_num > len(hash_key):
                c = exp(-pow(hash_key[point_num // len(hash_key)], 2) / 2)
            else:
                c = exp(-pow(hash_key[point_num-1], 2) / 2)
            gij = a * b * c

            if point_num >= len(hash_key):
                if hash_key[point_num // len(hash_key)] == 0:
                    edft = (hash_key[point_num // len(hash_key) + 2] * len(hash_key)) / gij
                else:
                    edft = (hash_key[point_num // len(hash_key)] * len(hash_key)) / gij
            else:
                if hash_key[point_num] == 0:
                    edft = (hash_key[point_num + 2] * len(hash_key)) / gij
                else:
                    edft = (hash_key[point_num] * len(hash_key)) / gij


            # 离散傅里叶变换
            # dft_ak.append(complex(X[feat_count][point_num], Y[feat_count][point_num]))
            #
            # # for l in range(len(X[feat_count])):
            # dft = np.fft.fft(dft_ak)

            ran_x += [X[feat_count][point_num] * gij]
            ran_y += [Y[feat_count][point_num] * gij]

            # Haar变换
            ap = ((ran_x[point_num] + ran_y[point_num]) / 2) * sqrt(2)
            dc = ((ran_x[point_num] - ran_y[point_num]) / 2) * sqrt(2)

            # 加密坐标值
            en_x += [ap * edft]
            en_y += [dc * edft]

            # 离散傅里叶反变换
            # idft_ak.append(complex(en_x[point_num],en_y[point_num]))
            # idft = np.fft.ifft(idft_ak)

            # Haar反变换
            iap = ((en_x[point_num] + en_y[point_num]) / 2) * sqrt(2)
            idc = ((en_x[point_num] - en_y[point_num]) / 2) * sqrt(2)
            A.append(iap), B.append(idc)


        En_X.append(A), En_Y.append(B)


    end_time = time()
    run_time = end_time - begin_time
    print('运行时间：', run_time)
    # print(En_X[0],En_X[-1])
    # print(En_Y[0],En_Y[-1])

    # print(En_X)
    # print(len(En_X))
    return En_X,En_Y



def write_encrytpion_shp(ori_shp,en_shp,user_key):

    En_X, En_Y = encryption_XY(ori_shp,user_key)

    ds = ogr.Open(ori_shp, 0)
    in_layer = ds.GetLayer(0)
    feature_num = in_layer.GetFeatureCount()

    driver = ogr.GetDriverByName('ESRI Shapefile')
    data_source = driver.CreateDataSource(en_shp)
    get_srs = in_layer.GetSpatialRef()

    out_layer = data_source.CreateLayer(en_shp, get_srs, ogr.wkbLineString)
    out_layer.CreateFields(in_layer.schema)
    out_defn = out_layer.GetLayerDefn()
    if os.access(en_shp, os.F_OK):
        driver.DeleteDataSource(en_shp)

    feat_count = 0
    for in_feature in in_layer:

        geom = in_feature.geometry()
        wkt = geom.ExportToWkt()

        pointstring = wkt[wkt.rfind('(') + 1:wkt.find(')')].split(',')
        pts_new = ''

        point_num = 0
        for point_value in pointstring:

            en_x = En_X[feat_count][point_num]
            en_y = En_Y[feat_count][point_num]

            pts_new += str(en_x) + ' ' + str(en_y) + ','

            point_num += 1

        en_wkt = wkt[0:wkt.find('(') + 1] + pts_new[:-1] + ')'
        out_feature = ogr.Feature(out_defn)
        en_Geometry = ogr.CreateGeometryFromWkt(en_wkt)

        out_feature.SetGeometry(en_Geometry)

        for i in range(1, in_feature.GetFieldCount()):
            value = in_feature.GetField(i)
            out_feature.SetField(i, value)

        out_layer.CreateFeature(out_feature)

        feat_count += 1
    data_source.Release()
    ds.Release()

def num_computer(ori_shp):
    ox,oy,feature_num = Get_XY_fromshp(ori_shp)
    vertice_num = 0
    for i in range(feature_num):
        vertice_num += len(ox[i])
    # print('初始x值：',ox)
    # print('初始y值：',oy)
    print('要素数量：',len(ox))
    print('顶点数量：',vertice_num)



if __name__ == '__main__':
    path = os.chdir(r'D:\2020\20200824 加密算法英文稿子\加密算法英文 稿子、数据与程序\test data')
    # Get_XY_fromshp('new_railways.shp')
    # encryption_XY('dunhuang_road.shp','cehuixueyuan2020')
    write_encrytpion_shp('gis_osm_railways_free_1.shp','en_gis_osm_railways_free_1.shp','cehuixueyuan2020')
    # num_computer('gis_osm_landuse_a_free_1_polygontoline.shp')
