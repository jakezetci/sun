# sun

# coordinates.py
Здесь сидит класс coordinates - который переводит сразу во все системы координат заданную точку. Точку можно задавать тремя способами - coordinates(x,y,z), coordinates(r, phi, theta, spherical=True), coordinates(r, lat, lon, latlon=True).

# field.py
Здесь сидят функции для модельных полей

# plots.py
Здесь большой файл для оформления графиков, пока оттуда нужна лишь функция sphere для построения сетки проекции сферы.

# lib.py

Здесь сидит большой класс grid, в котором отображается сетка, в которой можно задавать значение поля, а также функции для вычисления магнитного поля функцией Грина.

# testing_full.py

главный файл для тестирования - пока что он вычисляет поле на сетке [-30,30] [-30,30] по широте и долготе и строит график для этого поля
