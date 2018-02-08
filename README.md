# CellPi

__Under development, not ready for a use__

__English description will be uploaded with the first release__

This pipeline main idea is to take best-practice scRNA-seq data preprocessind tools and integrate them in a fully unsupervised pipeline "from GSE id to cluster's markers"

GSE_id->read counts will be based on the mixture of https://github.com/hms-dbmi/dropEst and https://github.com/CGATOxford/UMI-tools

Clustering step currently uses custom unsuperviced clustering. SC3 will be added as an alternative in a final version.

```
Инструкция:
1) Установить пакеты через installer.Rmd
2) Загрузить все необходимые функции из Counts2Exprs.Pmd
3) Загрузить данные ориентируясь на примеры из Load_data.Rmd
4) Запустить основной анализ по примерам из Simple_Analyse.md
5) Если функция seur_find_multiplication выдаёт NULL, значит разбить кластер на суб-кластеры не удалось

Аннотация исходной матрицы клеток должна иметь формат: <тип_клетки>.<номер_клетки> 
( "d7.122" "d0.484" "d2.39"  "d2.213")
Для расстановки номеров можно воспользоваться функцией types2names, 
пример использования которой можно найти в модуле Load_data.Rds в функции load_sce

Исходная аннотация клеток не используется при кластеризации но может быть полезна 
для оценки качества итогового разбиения

Переменная окружения "R_MAX_NUM_DLLS 255" может понадобиться при запуске в Windows, если возникает ошибка 
"превышение максимального числа загруженных библиотек"

```
