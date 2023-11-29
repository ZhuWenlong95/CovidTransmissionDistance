# CovidTransmissionDistance


## Linelist for the data and code
### Data
* data source<br>
This folder contains offline web pages with information on COVID-19 infections.<br>
Data on reported date and residential address of infections in COVID-19 outbreaks in cities of China were collected from city-level Health Commission/Center for Disease Control and Prevention.<br>
* variant info<br>
This folder contains offline web pages with information on each outbreak variant.<br>
Information on SARS-CoV-2 variants of each outbreak was collected from the reports of each city's Health Commission/Center for Disease Control and Prevention, research articles, or news.<br>
* CityCharacteristics.csv<br>
City characteristics used in the study.<br>
  * v.mean<br>
  Estimates of the mean transmission distance of each city.<br>
  * PopAll21 and CitySize<br>
  Number of permanent urban residents at the end of 2021 (PopAll21) and city size (CitySize).<br>
  From the national/city bureau of statistics and city’s seventh national population census bulletins.<br>
  * RoadLength<br>
  Total length of major urban roads of each city.<br>
  Calculated according to the map of roads in each city downloaded from Open Street Map (OSM, http://download.geofabrik.de/asia/china.html).<br>
  * BusStation<br>
  Total number of bus stations in each city.<br>
  From a bus tracing website (https://www.8684.cn).<br>
  * MobilityIndex<Br>
  The average of the intra-city travel intensity index from 1 March to 30 June 2023 in each city.<br>
  From the Baidu Migration Platform (https://qianxi.baidu.com).
### code
* TransmissionDistance 
