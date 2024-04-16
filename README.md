# Метод Шенкса (Shanks Transform) на C++

## Описание
Метод Шенкса, также известный как трансформация Ричардсона, представляет собой численный метод, применяемый для улучшения сходимости числовых последовательностей или рядов. Он находит применение в численном анализе и является методом улучшения точности численных вычислений.
В данном проекте на текущий момент реализован классический вариант метода Шенкса и Эпсилон алгоритм Винна.

Теория: [ссылка_на_теорию](https://drive.google.com/drive/folders/19KFEQhl9ZR4EE2zDFvi610bNdNBWfGIb?usp=sharing)

## Установка
В git bash

```
mkdir build
cd build
cmake ..
make
```

## Документация
Doxygen документация доступна [по ссылке](https://katerina-evdokimova.github.io/shanks-university/)


## Авторы
Большаков Михаил - Главный Программист (организация структуры кода, классический Шенкс) mike1024b@mail.ru

Евдокимова Екатерина - Программист (документация, комментарии, исправления) e_katerina.a_l@mail.ru

Пашков Борис - Программист (Эпсилон алгоритм, ряды), Помощь с Отчетом(Построение графиков, расчет времени выполнения) pashkovborya@gmail.com

Солониченко Злата - Теоретик, Отчет zlaaata_s@mail.ru

Кармалина Ольга - Теоретик, Отчет karmalinaolga@mail.ru

## Руководитель Проекта
Денис Васильевич Парфенов promasterden@yandex.ru

## Ссылки На Литературу
Статья Шенкса про его преобразования https://onlinelibrary.wiley.com/doi/abs/10.1002/sapm19553411

Публикации Винна (в частности есть статья про Эпсилон алгоритм Винна) https://mathresearch.utsa.edu/Legacy/Peter-Wynn/publications.html

Про Эпсилон алгоритм, примененный к монотонным и осцилирующим последовательностям https://www.sciencedirect.com/science/article/pii/S0377042700005616

Подробный анализ трансформации Шенкса и Эпсилон алгоритма https://www.researchgate.net/publication/327178717_The_genesis_and_early_developments_of_Aitken's_process_Shanks'_transformation_the_e-algorithm_and_related_fixed_point_methods
