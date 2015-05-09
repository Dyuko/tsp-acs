# tsp-acs

## About

This is an implementation of an Ant Colony System (ACS) to solve the Traveling Salesman Problem (TSP). It supports both sequential and parallel versions.

In ```doc/tsp-acs.pdf``` the algorithm is described in detail and a comparison of performance is discussed.

## Authors

* Felipe Ramirez

## Usage

```
./tsp-acs <data> <number_of_ants> <use_threads>
```
* ```data``` is an instance of the TSP
* ```number_of_ants``` is the number of ants to be routed
* ```use_threads``` whether ```0``` or ```1```