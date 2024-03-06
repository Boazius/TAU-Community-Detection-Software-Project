# Community Structure in Networks Detection

## Introduction
This project, developed in the C programming language for the Tel Aviv University Software Project Course in 2020/2021, implements an algorithm for detecting community structures, or clusters, in a network. These are groups of nodes within a network that have dense internal connections and sparser connections between groups. Identifying such structures is crucial in many fields, such as biology for understanding protein interactions.
This project provides a robust tool for community detection in complex networks, which is essential for the modular understanding of various real-world networks.

## Algorithm Description
The algorithm represents networks as graphs and uses modularity matrices to evaluate the quality of the divisions. It seeks to maximize modularity, which is the density of edges within groups versus between them.

### Division Algorithm
- The division algorithm starts by splitting the network into two communities, aiming for the maximum modularity.
- It then iteratively divides each community into smaller ones, recalculating the modularity matrix for each subdivision.
- The process includes optimization steps to refine the division and improve the modularity score.

### Leading Eigenpair
- The algorithm uses power iteration and matrix shifting techniques to find the leading eigenpair of the modularity matrix, crucial for the division process.

## Implementation Details
- The program, named `cluster`, takes two command-line arguments for input and output filenames.
- It reads a binary file describing the network graph and outputs a binary file with the community divisions.

## File Formats
- Both input and output files are binary, containing streams of integers representing nodes and edges for input and groups and their nodes for output.
