number_of_problems = 6;
number_of_runs = 10;

problem_name = {"Nonlinear", "Quad", "QuadFull", "HC", "m-RTAC", "Nagumo", "All"};

data = struct();

for i = 1:(number_of_problems)
  data.(problem_name{i}) = dlmread("output.csv", ',', [i-1, 1, i-1, number_of_runs]);
endfor

data.All = dlmread("output.csv", ',', [0, 1, number_of_problems, number_of_runs])

cdata = struct2cell(data);

stats = boxplot(cdata);

set(gca (), "xtick", 1:(number_of_problems+1), "xticklabel", problem_name);
xlabel('Benchmark Problems');
ylabel('parallelization (%), q');
title('Code parallelization estimate using benchmarking');

