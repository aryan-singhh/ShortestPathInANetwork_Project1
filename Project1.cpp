// C++ version of the given Java graph network program
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <queue>
#include <map>
#include <set>
#include <list>
#include <limits>
#include <algorithm>
using namespace std;

const double INFINITY_DIST = numeric_limits<double>::max();

class GraphException : public exception {
    string msg;
public:
    GraphException(string m) : msg(m) {}
    const char* what() const noexcept override {
        return msg.c_str();
    }
};

struct Vertex;

struct Edge {
    Vertex* dest;
    double distance;
    Edge(Vertex* d, double dist) : dest(d), distance(dist) {}
};

struct Vertex {
    string name;
    list<Edge> adj;
    Vertex* prev = nullptr;
    double dist = INFINITY_DIST;
    string status = "up";

    Vertex(string nm) : name(nm) {}
    void reset() {
        dist = INFINITY_DIST;
        prev = nullptr;
    }
};

struct Path {
    string str;
    double dist;
    Path(string s, double d) : str(s), dist(d) {}
    bool operator>(const Path& p) const {
        return dist > p.dist;
    }
};

struct Pair {
    string str1, str2;
    Pair(string s1, string s2) : str1(s1), str2(s2) {}
};

class Graph {
    map<string, Vertex*> vertexMap;
    list<Pair> edgeList;
    set<string> visited;

    Vertex* getVertex(const string& name) {
        if (vertexMap.find(name) == vertexMap.end())
            vertexMap[name] = new Vertex(name);
        return vertexMap[name];
    }

public:
    void addEdge(const string& src, const string& dest, double dist) {
        Vertex* v = getVertex(src);
        Vertex* w = getVertex(dest);
        auto it = find_if(v->adj.begin(), v->adj.end(), [&](Edge& e) { return e.dest->name == dest; });
        if (it != v->adj.end()) {
            it->distance = dist;
        } else {
            v->adj.emplace_back(w, dist);
        }
    }

    bool isEdgeDown(const string& src, const string& dest) {
        for (const auto& pair : edgeList)
            if (pair.str1 == src && pair.str2 == dest)
                return true;
        return false;
    }

    void print() {
        cout << "\n";
        set<string> nodes;
        for (const auto& p : vertexMap)
            nodes.insert(p.first);

        for (const string& node : nodes) {
            Vertex* v = vertexMap[node];
            cout << node;
            if (v->status != "up") cout << " down";
            cout << endl;

            map<string, double> neighbors;
            for (const Edge& e : v->adj)
                neighbors[e.dest->name] = e.distance;

            for (const auto& [dest, dist] : neighbors) {
                cout << "  " << dest << " " << dist;
                if (isEdgeDown(node, dest)) cout << " down";
                cout << endl;
            }
        }
    }

    void dijkstra(const string& start) {
        for (auto& [_, v] : vertexMap) v->reset();
        if (vertexMap[start]->status != "up") {
            cout << "Vertex is Down\n";
            return;
        }

        vertexMap[start]->dist = 0.0;
        priority_queue<Path, vector<Path>, greater<>> pq;
        pq.emplace(start, 0.0);

        while (!pq.empty()) {
            auto [uName, dist] = pq.top(); pq.pop();
            Vertex* u = vertexMap[uName];
            for (const Edge& e : u->adj) {
                Vertex* v = e.dest;
                if (v->status != "up" || isEdgeDown(u->name, v->name)) continue;
                double weight = e.distance;
                if (v->dist > u->dist + weight) {
                    v->dist = u->dist + weight;
                    v->prev = u;
                    pq.emplace(v->name, v->dist);
                }
            }
        }
    }

    void printPath(const string& dest) {
        if (vertexMap.find(dest) == vertexMap.end()) {
            cout << "Destination vertex not found\n";
            return;
        }
        Vertex* v = vertexMap[dest];
        if (v->dist == INFINITY_DIST) {
            cout << dest << " is unreachable\n";
            return;
        }

        vector<string> path;
        while (v) {
            path.push_back(v->name);
            v = v->prev;
        }
        reverse(path.begin(), path.end());
        for (const string& s : path) cout << s << " ";
        cout << fixed;
        cout.precision(2);
        cout << vertexMap[dest]->dist << endl;
    }

    void edgeDown(const string& src, const string& dest) {
        if (!isEdgeDown(src, dest))
            edgeList.emplace_back(src, dest);
    }

    void edgeUp(const string& src, const string& dest) {
        edgeList.remove_if([&](const Pair& p) { return p.str1 == src && p.str2 == dest; });
    }

    void vertexDown(const string& name) { vertexMap[name]->status = "down"; }
    void vertexUp(const string& name) { vertexMap[name]->status = "up"; }

    void deleteEdge(const string& src, const string& dest) {
        Vertex* v = getVertex(src);
        v->adj.remove_if([&](Edge& e) { return e.dest->name == dest; });
    }

    void reachableUtil(const string& node) {
        for (const Edge& e : vertexMap[node]->adj) {
            if (e.dest->status != "up" || isEdgeDown(node, e.dest->name)) continue;
            if (visited.insert(e.dest->name).second) {
                reachableUtil(e.dest->name);
            }
        }
    }

    void checkReachable() {
        set<string> nodes;
        for (auto& [k, v] : vertexMap)
            if (v->status == "up") nodes.insert(k);

        for (const string& node : nodes) {
            cout << node << endl;
            visited.clear();
            visited.insert(node);
            reachableUtil(node);
            for (const string& v : visited)
                if (v != node) cout << "  " << v << endl;
        }
    }

    void processCommand(const string& line) {
        stringstream ss(line);
        string cmd, arg1, arg2;
        double dist;
        ss >> cmd;
        if (cmd == "print") print();
        else if (cmd == "path") {
            ss >> arg1 >> arg2;
            dijkstra(arg1);
            printPath(arg2);
        }
        else if (cmd == "reachable") checkReachable();
        else if (cmd == "edgedown") { ss >> arg1 >> arg2; edgeDown(arg1, arg2); }
        else if (cmd == "edgeup") { ss >> arg1 >> arg2; edgeUp(arg1, arg2); }
        else if (cmd == "vertexdown") { ss >> arg1; vertexDown(arg1); }
        else if (cmd == "vertexup") { ss >> arg1; vertexUp(arg1); }
        else if (cmd == "deleteedge") { ss >> arg1 >> arg2; deleteEdge(arg1, arg2); }
        else if (cmd == "addedge") { ss >> arg1 >> arg2 >> dist; addEdge(arg1, arg2, dist); }
        else if (cmd == "quit") exit(0);
        else cout << "Not valid query\n";
    }

    void buildGraph(const string& filename) {
        ifstream infile(filename);
        string line;
        while (getline(infile, line)) {
            istringstream iss(line);
            string src, dest;
            double dist;
            if (!(iss >> src >> dest >> dist)) {
                cerr << "Skipping ill-formatted line: " << line << endl;
                continue;
            }
            addEdge(src, dest, dist);
            addEdge(dest, src, dist);
        }
    }
};

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <filename>\n";
        return 1;
    }

    Graph g;
    g.buildGraph(argv[1]);
    string line;
    while (getline(cin, line)) {
        g.processCommand(line);
    }
    return 0;
}
