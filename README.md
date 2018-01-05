data/ 目录下存放供程序的输入数据

doc/ 目录下存放电子版报告文档

src/ 目录下存放课程设计中编写的程序源代码

源程序由cmake管理,CMakeList.txt存放在src目录下

假设在 课程设计xinginx/build 目录下执行
cmake ../src
make -j4

编译后得到几个对应的可执行文件的基本使用为:

./EigenQuaternion

./ImageManipulate ../data/motto.png

./AccelAnalyse ../data/accel.txt

./GearAnalyse ../data/gear/*.txt

./FeatureExtraction ../data/1.jpg ../data/2.jpg

./PoseEstimation ../data/1.jpg ../data/2.jpg

./Triangulate ../data/1.jpg ../data/2.jpg
