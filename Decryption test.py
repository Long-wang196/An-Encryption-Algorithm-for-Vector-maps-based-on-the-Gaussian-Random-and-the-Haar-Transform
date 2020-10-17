import os,sys,ogr,gdal,hashlib,osr
from math import *
from time import *
import numpy as np

os.chdir(r'D:\2020\20200824 加密算法英文稿子\加密算法英文 稿子、数据与程序\test data')
def gene_hashkey(user_key):
    # 密钥K的生成
    sha = hashlib.sha512()
    sha.update(user_key.encode('utf-8'))
    hash_key = sha.hexdigest()
    value_list = [int(i, 16) for i in list(hash_key)]
    return value_list

def Get_enXY_fromshp(en_shp):

    X, Y = [], []

    ds = ogr.Open(en_shp,0)
    in_layer = ds.GetLayer(0)
    feature_num = in_layer.GetFeatureCount()

    for in_feature in in_layer:

        geom = in_feature.geometry()
        wkt = geom.ExportToWkt()

        pointstring = wkt[wkt.find('(') + 1:wkt.rfind(')')].split(',')

        x, y = [],[]
        for point_value in pointstring:

            if point_value == '':
                continue
            else:
                xy = point_value.split(' ')

                x += [float(xy[0])]
                y += [float(xy[1])]

        X.append(x), Y.append(y)
    # print(X)
    # print(len(X),feature_num)
    del ds
    return X,Y,feature_num


def decryption_XY(en_shp,user_key):
    hash_key = gene_hashkey(user_key)
    X, Y, feature_num = Get_enXY_fromshp(en_shp)
    DEn_X, DEn_Y = [], []

    for feat_count in range(feature_num):

        dft_ak, idft_ak = [], []
        al, alx, aly = [], [], []
        den_x, den_y = [], []
        ran_x, ran_y = [], []
        idft_Al, idft_Alx, idft_Aly = [], [], []
        for point_num in range(len(X[feat_count])):
            a = ((feat_count + 1) * (point_num + 1)) / feature_num
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
            # dft_ak += [complex(X[feat_count][point_num],Y[feat_count][point_num])]
            # dft = np.fft.fft(dft_ak)

            # Haar变换
            ap = ((X[feat_count][point_num] + Y[feat_count][point_num]) / 2) * sqrt(2)
            dc = ((X[feat_count][point_num] - Y[feat_count][point_num]) / 2) * sqrt(2)

            # 解密坐标值
            den_x += [ap / edft]
            den_y += [dc / edft]

            # 离散傅里叶逆变换
            # idft_ak += [complex(den_x[point_num],den_y[point_num])]
            # idft = np.fft.ifft(idft_ak)

            # Haar反变换
            iap = ((den_x[point_num] + den_y[point_num]) / 2) * sqrt(2)
            idc = ((den_x[point_num] - den_y[point_num]) / 2) * sqrt(2)

            ran_x += [iap / gij]
            ran_y += [idc / gij]
            # alx.append(ran_x), aly.append(ran_y)

        DEn_X.append(ran_x), DEn_Y.append(ran_y)

    # print(DEn_X)
    # print(DEn_Y)
    # print(len(DEn_X),len(DEn_Y))
    return DEn_X,DEn_Y

def write_decryption_shp(en_shp,de_shp,user_key):
    De_X, De_Y = decryption_XY(en_shp, user_key)

    ds = ogr.Open(en_shp)
    in_layer = ds.GetLayer()
    feature_num = in_layer.GetFeatureCount()
    get_srs = in_layer.GetSpatialRef()

    driver = ogr.GetDriverByName('ESRI Shapefile')
    data_source = driver.CreateDataSource(de_shp)
    out_layer = data_source.CreateLayer(de_shp, get_srs, ogr.wkbLineString)
    out_layer.CreateFields(in_layer.schema)
    out_defn = out_layer.GetLayerDefn()
    if os.access(de_shp, os.F_OK):
        driver.DeleteDataSource(de_shp)

    feat_count = 0
    # de_wkt = []
    for in_feature in in_layer:

        geom = in_feature.geometry()
        wkt = geom.ExportToWkt()

        pointstring = wkt[wkt.find('(') + 1:wkt.rfind(')')].split(',')
        pts_new = ''

        point_num = 0
        for point_value in pointstring:
            de_x = De_X[feat_count][point_num]
            de_y = De_Y[feat_count][point_num]

            pts_new += str(de_x) + ' ' + str(de_y) + ','

            point_num += 1

        de_wkt = wkt[0:wkt.find('(') + 1] + pts_new[:-1] + ')'
        out_feature = ogr.Feature(out_defn)
        de_Geometry = ogr.CreateGeometryFromWkt(de_wkt)

        out_feature.SetGeometry(de_Geometry)

        for i in range(1, in_feature.GetFieldCount()):
            value = in_feature.GetField(i)
            out_feature.SetField(i, value)

        out_layer.CreateFeature(out_feature)
    #
        feat_count += 1

    del data_source
    # print(de_wkt)
    # print(in_layer.GetSpatialRef())

if __name__ == '__main__':
    path = os.chdir(r'D:\2020\20200824 加密算法英文稿子\加密算法英文 稿子、数据与程序\test data')
    # Get_enXY_fromshp('en_road.shp')
    # decryption_XY('en_road.shp','cehuixueyuan2020')
    write_decryption_shp('en_gis_osm_railways_free_1.shp','de_gis_osm_railways_free_1.shp','cehuixueyuan2020')