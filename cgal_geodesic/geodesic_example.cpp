#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh_shortest_path.h>

#include <iostream>
#include <fstream>
#include <exception>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     Kernel;
typedef Kernel::Point_3                                         Point_3;
typedef CGAL::Surface_mesh<Point_3>                             Mesh;
typedef CGAL::Surface_mesh_shortest_path_traits<Kernel, Mesh>  Traits;
typedef CGAL::Surface_mesh_shortest_path<Traits>                Shortest_path;

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " input.off source_vertex_index" << std::endl;
        return 1;
    }

    const char* filename = argv[1];
    int source_idx = std::stoi(argv[2]);

    Mesh mesh;
    std::ifstream input(argv[1]);
    if (!input) {
        std::cerr << "Error: Could not open file." << std::endl;
        return 1;
    }
    if (!(input >> mesh)) {
        std::cerr << "Error: Could not parse mesh from file." << std::endl;
        return 1;
    }
    if (mesh.is_empty()) {
        std::cerr << "Error: Mesh is empty after reading." << std::endl;
        return 1;
    }

    if (source_idx < 0 || source_idx >= static_cast<int>(num_vertices(mesh))) {
        std::cerr << "Error: source vertex index out of range." << std::endl;
        return 1;
    }

    Mesh::Vertex_index source_v(source_idx);

    Shortest_path shortest_paths(mesh);
    shortest_paths.add_source_point(source_v);

    // Compute and print distance to all other vertices
    std::ofstream out("distances.txt");

    for (auto vd : mesh.vertices()) {
    double d = shortest_paths.shortest_distance_to_source_points(vd).first;
    out << vd.idx() << " " << d << "\n";
    }

    out.close();

    // HERE STARTS THE PATH VISUALIZATION

    // Choose a target vertex for visualization of path (not the source itself)
    int target_idx = 2302;
    if (target_idx == source_idx) {
        target_idx = 1; // Just pick the next vertex if equal
    }
    if (target_idx < 0 || target_idx >= static_cast<int>(num_vertices(mesh))) {
        std::cerr << "Error: target vertex index out of range." << std::endl;
        return 1;
    }

    Mesh::Vertex_index target_v(target_idx);

    // Extract shortest path points from target back to source
    std::vector<Point_3> path_points;
    shortest_paths.shortest_path_points_to_source_points(target_v, std::back_inserter(path_points));

    // Write path points to plain text file
    std::ofstream path_out("path.txt");
    if (!path_out) {
        std::cerr << "Error: cannot open path.txt for writing." << std::endl;
        return 1;
    }

    for (const auto& p : path_points) {
        path_out << p.x() << " " << p.y() << " " << p.z() << "\n";
    }

    path_out.close();
    std::cout << "Shortest path written to path.txt (from vertex " << target_idx
              << " to source vertex " << source_idx << ")" << std::endl;

    return 0;
}