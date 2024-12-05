#include "GpuSpec.h"
#include <gtest/gtest.h>

TEST(GpuSpecTest, PrintTest) {
    GpuSpec gpu("NVIDIA", 8192);
    testing::internal::CaptureStdout();
    gpu.Print();
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_EQ(output, "GPU Model: NVIDIA, Memory: 8192 MB\n");
}

TEST(GpuSpecTest, ImportExportTest) {
    GpuSpec gpu("AMD", 4096);
    gpu.Export("gpu_test.txt");
    GpuSpec importedGpu;
    importedGpu.Import("gpu_test.txt");
    testing::internal::CaptureStdout();
    importedGpu.Print();
    std::string output = testing::internal::GetCapturedStdout();
    EXPECT_EQ(output, "GPU Model: AMD, Memory: 4096 MB\n");
}
