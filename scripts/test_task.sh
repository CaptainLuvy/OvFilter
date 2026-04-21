#!/bin/bash

echo ">>> [Step 1] 任务初始化..."
sleep 2

echo ">>> [Step 2] 正在处理数据 (模拟耗时)..."
for i in {1..3}; do
    echo "    处理进度: $i/3..."
    sleep 1
done

# 如果参数是 "fail"，则模拟报错
if [ "$1" == "fail" ]; then
    echo "!!! [Error] 发生严重错误！"
    echo "!!! [Error] 模拟脚本异常退出。"
    exit 1
fi

echo ">>> [Step 3] 保存结果..."
sleep 1
echo ">>> [Success] 所有步骤已完成！"
exit 0
