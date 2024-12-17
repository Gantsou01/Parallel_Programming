#include "GpuSpec.h"
#include "CpuSpec.h"
#include "RpmSpec.h"
#include "LanSpec.h"
#include "ClusterNode.h"
#include "Cluster.h"

int main() {
    GpuSpec gpu("NVIDIA", 8192);
    CpuSpec cpu("Intel Xeon", 16);
    RpmSpec rpm("SSD", 512);
    LanSpec lan("Ethernet", 10000);

    ClusterNode node(gpu, cpu, rpm, lan);
    node.Print();

    Cluster cluster;
    cluster.AddNode(node);
    cluster.Print();

    // Экспорт и импорт данных
    cluster.Export("cluster_data.txt");
    cluster.Import("cluster_data.txt");

    return 0;
}
