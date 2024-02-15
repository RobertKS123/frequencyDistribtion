using System;
using System.Collections.Generic;
using System.Drawing;
using static System.Runtime.InteropServices.JavaScript.JSType;

namespace frequency_distribution
{
    class Program
    {
        static void Main(string[] args)
        {
            GetColors getColors = new GetColors();
            getColors.getTowers();
            getColors.calculateColors();
            //getColors.calculateColorsGreedy();
        }
    }

    class GetColors
    {
        private int n;
        private List<Tower> towers;
        private double[][] distanceMatrix;
        private double maxDistance = 0;
        private List<int> bestSolutionColors;
        private double bestSolution = double.MaxValue;
        private bool[][] adjacencyMatrix;

        // user enters towers or use default provided towers
        public void getTowers(bool defaultTowers = true)
        {
            if (defaultTowers)
            {
                towers = new List<Tower>
                {
                    new Tower("A", 536660, 183800, -0.03098, 51.53657),
                    new Tower("B", 537032, 184006, -0.02554, 51.53833),
                    new Tower("C", 537109, 183884, -0.02448, 51.53721),
                    new Tower("D", 537110, 184695, -0.02415, 51.5445),
                    new Tower("E", 537206, 184685, -0.02277, 51.54439),
                    new Tower("F", 537248, 185016, -0.02204, 51.54735),
                    new Tower("G", 537250, 185020, -0.02201, 51.54739),
                    new Tower("H", 537267, 184783, -0.02185, 51.54525),
                    new Tower("I", 537269, 183451, -0.02234, 51.53328),
                    new Tower("J", 537270, 184140, -0.02206, 51.53948),
                    new Tower("K", 537356, 184927, -0.02052, 51.54653),
                    new Tower("L", 537380, 184727, -0.02025, 51.54472),
                    new Tower("M", 537458, 184495, -0.01921, 51.54262),
                    new Tower("N", 537604, 184134, -0.01725, 51.53934),
                    new Tower("O", 537720, 184057, -0.01561, 51.53862),
                    new Tower("P", 537905, 184591, -0.01273, 51.54337),
                    new Tower("Q", 537910, 184441, -0.01272, 51.54202),
                    new Tower("R", 537953, 184295, -0.01216, 51.5407),
                    new Tower("S", 538050, 184245, -0.01078, 51.54023)
                };
            }
            else
            {
                towers = new List<Tower>();

                Console.WriteLine("Enter tower data (name, easting, northing, long, lat), or -1 to finish:");

                while (true)
                {
                    string input = Console.ReadLine();

                    if (input == "-1")
                        break;

                    string[] parts = input.Split(',');

                    if (parts.Length != 5)
                    {
                        Console.WriteLine("Invalid input format. Please try again.");
                        continue;
                    }

                    string name = parts[0].Trim();
                    double easting, northing, longitude, latitude;

                    if (!double.TryParse(parts[1].Trim(), out easting) ||
                        !double.TryParse(parts[2].Trim(), out northing) ||
                        !double.TryParse(parts[3].Trim(), out longitude) ||
                        !double.TryParse(parts[4].Trim(), out latitude))
                    {
                        Console.WriteLine("Invalid numeric format. Please try again.");
                        continue;
                    }

                    towers.Add(new Tower(name, easting, northing, longitude, latitude));
                }

                Console.WriteLine("Getting Frequency distribution");
            }
        }

        // All nessisary steps to calcutale tower frequencies 
        public void calculateColors()
        {

            n = towers.Count;

            // Create a distance matrix so distance calculations only need to happen once
            getDistances();

            // get a distance slightly bigger than the largest distance between two points for a negative reward
            maxDistance = distanceMatrix.SelectMany(subArray => subArray ?? Enumerable.Empty<double>()).Max() + 1;

            // List to store the best solution
            bestSolutionColors = new List<int>();

            // List for the current itteration in solution space
            List<int> currentColors = new List<int>();

            // Find optimal solution
            assignFrequencies(0, 0, currentColors, 0);

            // Output best solution
            Console.WriteLine("Best solution colors:");
            for (int i = 0; i < n; i++)
            {
                Console.Write(towers[i].Name + (int)(bestSolutionColors[i] + 110) + " ");
            }
        }

        // Recursivly expolre the solution space to find best solution in a DFS manner
        // Loops throught the Cell Towers sequentially
        // Give each a color and check if better than current best solution
        private void assignFrequencies(int i, double value, List<int> currentColors, int startColor)
        {
            //base case, when we hit a leaf node
            if (i >= n)
            {
                // is this the best leaf update the best potential solution
                if (value < bestSolution)
                {
                    bestSolution = value;
                    bestSolutionColors = new List<int>(currentColors);
                }
                return;
            }

            Tower t = towers[i];
            // starts loop at previous frequency +1
            for (int color = 0; color < 6; color++)
            {
                int current_color = (color + startColor) % 6;

                // only add to sum if colors are the same
                // by the bottom if the tree should have the total 
                double newValue = value;

                // get value of current state 
                foreach (int c in currentColors)
                {
                    if (current_color == c)
                    {
                        // the bigger the distance the better 
                        newValue += maxDistance - distanceMatrix[color][c];
                    }
                }
                if (newValue > bestSolution)
                {
                    // we already have a guarantied better solution so continue
                    continue;
                }

                currentColors.Add(current_color);
                assignFrequencies(i + 1, newValue, currentColors, current_color + 1);
                //backtracks by removing tower t from list and re-adding it with a different color
                currentColors.RemoveAt(currentColors.Count - 1);
            }
        }

        // get the distance between each tower 
        // matrix is symmetrical 
        private void getDistances()
        {
            // initalize the matrix
            distanceMatrix = new double[n][];
            for (int i = 0; i < n; i++)
            {
                distanceMatrix[i] = new double[n];
            }
            // do calcutalions
            for (int i = 0; i < n; i++)
            {
                for (int j = i + 1; j < n; j++)
                {
                    double d = getDistance(towers[i], towers[j]);
                    distanceMatrix[i][j] = d;
                    distanceMatrix[j][i] = d;
                }
            }
        }

        // Calculate the distance using the Haversine method to account for latitude and longitude
        private double getDistance(Tower a, Tower b)
        {
            const double earthRadius = 6371;

            var lat1 = degreesToRadians(a.Lat);
            var long1 = degreesToRadians(a.Long);
            var lat2 = degreesToRadians(b.Lat);
            var long2 = degreesToRadians(b.Long);

            var dLat = lat2 - lat1;
            var dLond = long2 - long1;

            var haverA = Math.Sin(dLat / 2) * Math.Sin(dLat / 2) +
                       Math.Cos(lat1) * Math.Cos(lat2) *
                       Math.Sin(dLond / 2) * Math.Sin(dLond / 2);
            var haverC = 2 * Math.Atan2(Math.Sqrt(haverA), Math.Sqrt(1 - haverA));

            return earthRadius * haverC;
        }

        // converts degrees to radians
        private double degreesToRadians(double d)
        {
            return d * Math.PI / 180;
        }

        // Greedy appoch to solution using a coloring alogrithm
        public void calculateColorsGreedy()
        {
            n = towers.Count;

            // Create a distance matrix so distance calculations only need to happen once
            getDistances();

            // get a distance slightly bigger than the largest distance between two points for a negative reward
            maxDistance = distanceMatrix.SelectMany(subArray => subArray ?? Enumerable.Empty<double>()).Max() + 1;

            // get the adjacency Matrix
            AddAdjacencies(0.5); // change to a better value

            greedyFrequencies();

            foreach (Tower t in towers)
            {
                Console.WriteLine(t.Name + " " + t.Long + " " + t.Lat + " " + (int)(t.Frequency +110));
            }
        }

        // Creates a graph using an adjasencies based on a proximity threshhold
        // towers need at least one neighbor 
        // towers can have a maximum of 5 neighbors (5 closest are choosen)
        void AddAdjacencies(double threshold)
        {
            foreach (Tower tower in towers)
            {
                FindClosestTowers(tower, threshold);
            }
        }

        void FindClosestTowers(Tower currentTower, double threshold)
        {
            List<Tower> closestTowers = new List<Tower>();
            Tower closestTower = null;
            double closestDistance = double.MaxValue;

            foreach (var tower in towers)
            {
                if (tower == currentTower)
                    continue;

                double distance = getDistance(currentTower,tower);

                if (distance <= threshold)
                {
                    closestTowers.Add(tower);
                }
                if (distance < closestDistance)
                {
                    closestDistance = distance;
                    closestTower = tower;
                }
            }
            if (closestTowers.Count > 0)
            {
                closestTowers.Sort((t1, t2) => getDistance(currentTower, t1).CompareTo(getDistance(currentTower, t2)));
                for (int i = 0; i < Math.Min(closestTowers.Count , 5); i++)
                {
                    currentTower.AddNeighbor(closestTowers[i]);
                }
            } 
            else
            {
                currentTower.AddNeighbor(closestTower);
            }
        }

        // Assign each point a frequency in a greedy fashon
        void greedyFrequencies()
        {

            // start with a tower with max neighbors
            Tower start = towers.OrderByDescending(t => t.Neighbors.Count).First();
            start.Frequency = 0;

            // keep track of used frequencies 
            int[] colors = new int[6];

            foreach (var neighbor in start.Neighbors)
            {
                // Assign a color different from the start tower
                neighbor.Frequency = (start.Frequency + 1) % 6; 
                colors[neighbor.Frequency] = 1;
            }

            // for each tower except the start tower 
            foreach (var tower in towers.Where(t => t != start))
            {
                // fidn the available colors for this tower
                bool[] availableColors = new bool[6];
                foreach (var neighbor in tower.Neighbors)
                {
                    // if neighbor has a color mark color as used
                    if (neighbor.Frequency != -1)
                    {
                        availableColors[neighbor.Frequency] = true;
                    }
                }

                // Find the first available color
                int color = Array.IndexOf(availableColors, false);
                if (color == -1)
                {
                    // If all colors are used, find the first unused color
                    color = Array.IndexOf(colors, 0);
                    colors[color] = 1; // Mark the color as used
                }

                // Assign the frequency to the tower
                tower.Frequency = color;
            }
        }
    }

    class Tower
    {
        public string Name { get; }
        public double Easting { get; }
        public double Northing { get; }
        public double Long { get; }
        public double Lat { get; }

        public List<Tower> Neighbors = new List<Tower>();
        public int Frequency = -1; // color

        public Tower(string name, double easting, double northing, double longitude, double lat)
        {
            Name = name;
            Easting = easting;
            Northing = northing;
            Long = longitude;
            Lat = lat;
        }

        public void AddNeighbor(Tower neighbor)
        {
            Neighbors.Add(neighbor);
        }
    }

}
