# Block-Gauss
Этот проект - мое решение приведенной ниже задачи.

Статус проекта: в разработке


## Содержание
- [Задача](#задача)
- [Использование](#использование)
- [Тестирование](#тестирование)

## Задача
- Требуется написать блочный алгоритм нахождения решения линейной системы уравнений методом Гаусса с выбором главного элемента по всей матрице.
Инициализация матрицы происходит либо по формулам, либо из файла. Инициализация вектор-столбца В происходит после инициализации матрицы А путем суммирования ее нечетных столбцов.
В конце программа выдает две невязки и время работы двух функций: функции, вычисляющей невязку и ищущую решение СЛУ.

## Использование
Программа работает в двух режимах: считывание из файла(1) и инициализация по формуле(2). После сборки программы утилитой make:
- (1)
```sh
$ ./a.out n m r 0 filename.txt
```
- (2)
```sh
$ ./a.out n m r s
```
Где n - размер матрицы, m - размер блока, r - максимальное количество отображаемых элементов, s - номер формулы, filename.txt - название файла, где лежит матрица A.

## Тестирование

В директории tests находится 2 shell скрипта. Чтобы запустить тесты:
```sh
$ bash filename.sh
```