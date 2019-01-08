#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <vector>

// struct Point {
//   float x{0}, y{0};
//   vector<float> x;
// };

using Point = std::vector<float>;
using DataFrame = std::vector<Point>;

float square(float value) {
  return value * value;
}

float squared_l2_distance(Point first, Point second) {
  int inc=0;
  float sum = 0.0;
  for (inc = 0 ; inc < (int)first.size() ; inc++)
    {
      sum = sum + square(first[inc] - second[inc]);
    }
  return sum;
}
std::vector<size_t> k_means_classif(const DataFrame& data, const DataFrame& mean)
{
    std::vector<size_t> to_return;
    for (size_t point = 0; point < data.size(); point++) 
    {
        float best_distance = 1000.0;
        size_t best_cluster = 0;
        for (size_t cluster = 0; cluster < mean.size(); cluster++) 
        {
            const float distance =
                squared_l2_distance(data[point], mean[cluster]);
            if (distance < best_distance) 
            {
                best_distance = distance;
                best_cluster = cluster;
            }
        }
        to_return.push_back(best_cluster);
    } 
    return to_return;
}

DataFrame k_means(const DataFrame& data,
                  size_t k,
                  size_t number_of_iterations) {
  static std::random_device seed;
  static std::mt19937 random_number_generator(seed());
  std::uniform_int_distribution<size_t> indices(0, data.size() - 1);

  // Pick centroids as random points from the dataset.
  DataFrame means(k);
  for (auto& cluster : means) {
    cluster = data[indices(random_number_generator)];
  }

  std::vector<size_t> assignments(data.size());
  for (size_t iteration = 0; iteration < number_of_iterations; ++iteration) {
    std::cerr << "." ;
    // Find assignments.
    for (size_t point = 0; point < data.size(); ++point) {
      auto best_distance = std::numeric_limits<float>::max();
      size_t best_cluster = 0;
      for (size_t cluster = 0; cluster < k; ++cluster) {
        const float distance =
            squared_l2_distance(data[point], means[cluster]);
        if (distance < best_distance) {
          best_distance = distance;
          best_cluster = cluster;
        }
      }
      assignments[point] = best_cluster;
    }

    // Sum up and count points for each cluster.
    DataFrame new_means(k);
    std::vector<size_t> counts(k, 0);
    for (size_t point = 0; point < data.size(); ++point) {
      const auto cluster = assignments[point];
      for (int inc = 0; inc < (int)new_means[cluster].size(); inc++)
      {
          new_means[cluster].at(inc) += data[point].at(inc);
//       new_means[cluster].y += data[point].y;
      }
      counts[cluster] += 1;
    }

    // Divide sums by counts to get new centroids.
    for (size_t cluster = 0; cluster < k; ++cluster) {
      // Turn 0/0 into 0/1 to avoid zero division.
      const auto count = std::max<size_t>(1, counts[cluster]);
      for (int inc = 0; inc < (int)new_means[cluster].size(); inc++)
      {
          means[cluster].at(inc) = new_means[cluster].at(inc) / count;
//       new_means[cluster].y += data[point].y;
      }
//       means[cluster].x = new_means[cluster].x / count;
//       means[cluster].y = new_means[cluster].y / count;
    }
  }

    std::cerr << std::endl;
  return means;
}

int main(int argc, const char* argv[]) {
  if (argc < 3) {
    std::cerr << "usage: k_means <data-file> <cluster-output-file> <data-clustered-output-file> <k> <vector-size> [iterations] [runs]"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }

  const auto k = std::atoi(argv[4]);
  const auto size_vector = std::atoi(argv[5]);
  const auto iterations = (argc >= 7) ? std::atoi(argv[6]) : 300;
  const auto number_of_runs = (argc >= 8) ? std::atoi(argv[7]) : 10;
  DataFrame data;
  std::ifstream stream(argv[1]);
  if (!stream) {
    std::cerr << "Could not open file: " << argv[1] << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::string line;
  long l_input_cpt=0;
  while (std::getline(stream, line)) {
    Point point;
    std::istringstream line_stream(line);
    size_t label;
    l_input_cpt++;
    if (l_input_cpt % 10000 == 0) std::cerr << ".";
    if (l_input_cpt % 100000 == 0) std::cerr << "(" << l_input_cpt << ")";

    for (int inc = 0; inc < size_vector ; inc++)
    {
      float datapoint=0.0;
      line_stream >> datapoint;
      point.push_back(datapoint);
    }
    data.push_back(point);
  }
  std::cerr << std::endl;

  DataFrame means;
  double total_elapsed = 0;
  for (int run = 0; run < number_of_runs; ++run) {
    std::cerr << "Run number: " << run << std::endl;
    const auto start = std::chrono::high_resolution_clock::now();
    means = k_means(data, k, iterations);
    const auto end = std::chrono::high_resolution_clock::now();
    const auto duration =
        std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    total_elapsed += duration.count();
  }
  std::cerr << "Took: " << total_elapsed / number_of_runs << "s ("
            << number_of_runs << " runs)" << std::endl;

  std::ofstream ocstream(argv[2]);
  if (!ocstream) {
    std::cerr << "Could not open file: " << argv[2] << std::endl;
    std::exit(EXIT_FAILURE);
  }
  
  for (auto& mean : means) 
  {
      for (int inc = 0; inc < (int)mean.size(); inc++)
      {
          ocstream << mean.at(inc) << " ";
      }
      ocstream << std::endl;
  }
  ocstream.close();
  
  std::vector <size_t> classifs = k_means_classif(data,means);
  std::ofstream odstream(argv[3]);
  if (!odstream) {
    std::cerr << "Could not open file: " << argv[3] << std::endl;
    std::exit(EXIT_FAILURE);
  }
  for (int inc = 0; inc < (int)classifs.size(); inc++)
  {
      odstream << classifs.at(inc) << std::endl;
  }
  odstream.close();
}
