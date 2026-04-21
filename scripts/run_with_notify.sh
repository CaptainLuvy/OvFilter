#!/bin/bash
# ==============================================================================
# 脚本名称: run_with_notify.sh
# 功能描述: 
#   1. 包装器：执行任意命令，并在命令结束后发送通知（包含耗时、状态、最后几行输出）。
#   2. 通知发送：内置支持 PushPlus/微信, 钉钉, Bark/iOS (官方或自建服务端)。
#      (v4.0: 终端输出全英文避免乱码，通知内容由 Python 生成确保中文正确)
#
# 用法: 
#   bash run_with_notify.sh [任务名称] <实际要运行的命令...>
# ==============================================================================

# --- 1. 配置区域 ---

# [PushPlus Token]
PUSHPLUS_TOKEN=""

# [Bark URL]
# ? 请务必删除下面的示例 Key，换成您自己的！
# 支持官方服务器 (https://api.day.app/...) 或自建服务器 (http://IP:Port/...)
BARK_URL="https://api.day.app/E2qf5FeXvWx6CJp2puN9Qm/"

# [钉钉 Webhook]
DINGTALK_WEBHOOK=""

# ------------------------------------------------

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 [Task_Name] <Command...>"
    echo "   or: $0 <Command...>"
    exit 1
fi

# 智能判断:
# 1. 如果 $1 是已知的命令(在PATH中)、文件、路径，或者常见的解释器(bash/python)，则认为没有提供任务名。
# 2. 否则，认为 $1 是任务名。

FIRST_ARG="$1"
CMD_ARRAY=()
TASK_NAME=""

# 检查 $1 是否是命令/文件/解释器
IS_CMD=false
if [[ "$FIRST_ARG" == "bash" ]] || [[ "$FIRST_ARG" == "sh" ]] || [[ "$FIRST_ARG" == "python"* ]] || \
   [[ -x "$FIRST_ARG" ]] || [[ -f "$FIRST_ARG" ]] || \
   [[ "$FIRST_ARG" == ./* ]] || [[ "$FIRST_ARG" == /* ]] || \
   command -v "$FIRST_ARG" &> /dev/null; then
    IS_CMD=true
fi

if [ "$IS_CMD" = true ]; then
    # 情况1: 用户直接输入了命令，没有提供任务名
    # 自动推断任务名
    if [[ "$FIRST_ARG" == "bash" ]] || [[ "$FIRST_ARG" == "sh" ]] || [[ "$FIRST_ARG" == "python"* ]]; then
        # 如果是 bash script.sh，取 script.sh 为名
        if [ -n "$2" ]; then
            TASK_NAME=$(basename "$2")
        else
            TASK_NAME=$(basename "$FIRST_ARG")
        fi
    else
        TASK_NAME=$(basename "$FIRST_ARG")
    fi
    CMD_ARRAY=("$@")
else
    # 情况2: 用户提供了任务名
    TASK_NAME="$1"
    shift
    CMD_ARRAY=("$@")
fi

# 智能添加解释器 (针对 .sh 和 .py 文件)
# 检查实际执行的命令部分
REAL_CMD="${CMD_ARRAY[0]}"

if [[ "$REAL_CMD" == *".sh" ]] && [[ "$REAL_CMD" != "bash" ]] && [[ "$REAL_CMD" != "sh" ]]; then
     CMD_ARRAY=("bash" "${CMD_ARRAY[@]}")
elif [[ "$REAL_CMD" == *".py" ]] && [[ "$REAL_CMD" != "python"* ]]; then
     if command -v python3 &> /dev/null; then
         CMD_ARRAY=("python3" "${CMD_ARRAY[@]}")
     else
         CMD_ARRAY=("python" "${CMD_ARRAY[@]}")
     fi
fi

# 生成用于显示的命令字符串
CMD_STR="${CMD_ARRAY[*]}"

LOG_TMP=$(mktemp)

echo "=========================================================="
echo "   [Monitor] Task Started: $TASK_NAME"
echo "   [Command] $CMD_STR"
echo "=========================================================="

START_TIME=$(date +%s)

# 执行命令
set +e
"${CMD_ARRAY[@]}" 2>&1 | tee "$LOG_TMP"
EXIT_CODE=${PIPESTATUS[0]}
set -e

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

# 计算耗时 (分:秒)
MIN=$((DURATION / 60))
SEC=$((DURATION % 60))
TIME_STR="${MIN}m ${SEC}s"

echo "=========================================================="
echo "   [Monitor] Task Finished (Duration: $TIME_STR)"
echo "   [Notify] Sending notification..."

# 优先使用 python3
PY_CMD="python3"
if ! command -v python3 &> /dev/null; then
    PY_CMD="python"
fi

# 生成临时 Python 脚本 (使用 Unicode 转义避免任何编码问题)
NOTIFY_SCRIPT="notify_temp_$$.py"

cat <<EOF > "$NOTIFY_SCRIPT"
# -*- coding: utf-8 -*-
import urllib.request
import urllib.parse
import json
import sys
import os
import codecs

# 强制 stdout/stderr 使用 utf-8 (解决 Windows 控制台乱码)
sys.stdout = codecs.getwriter("utf-8")(sys.stdout.detach())
sys.stderr = codecs.getwriter("utf-8")(sys.stderr.detach())

# 获取参数
if len(sys.argv) < 6:
    print("Error: Missing arguments")
    sys.exit(1)

exit_code = int(sys.argv[1])
task_name = sys.argv[2]     # 可能包含中文
duration_str = sys.argv[3]
cmd_str = sys.argv[4]       # 可能包含中文
log_file = sys.argv[5]

pushplus_token = '$PUSHPLUS_TOKEN'
bark_url = '$BARK_URL'
dingtalk_url = '$DINGTALK_WEBHOOK'

# 读取日志摘要
log_tail = ""
try:
    if os.path.exists(log_file):
        # 1. 优先尝试 UTF-8 (最通用)
        try:
            with open(log_file, 'r', encoding='utf-8') as f:
                lines = f.readlines()
        except UnicodeDecodeError:
            # 2. 如果失败，尝试 GBK (兼容 Windows 中文环境)
            with open(log_file, 'r', encoding='gbk', errors='replace') as f:
                lines = f.readlines()
        
        log_tail = "".join(lines[-15:])
except Exception as e:
    log_tail = f"(Log read error: {e})"

# 构造通知内容 (使用 Unicode 转义确保中文正确)
# \u2705 = ?, \u274c = ?, \u23f1 = ?, \u6210\u529f = 成功, \u5931\u8d25 = 失败
# \u8017\u65f6 = 耗时, \u547d\u4ee4 = 命令, \u6700\u540e\u8f93\u51fa = 最后输出

if exit_code == 0:
    title = f"[\u2705 \u6210\u529f] {task_name}"
    status_text = "Success"
    status_icon = "\u2705"
else:
    title = f"[\u274c \u5931\u8d25] {task_name}"
    status_text = f"Failed (Code {exit_code})"
    status_icon = "\u274c"

msg = f"""{status_icon} Status: {status_text}
\u23f1 Duration: {duration_str}
\uD83D\uDCDD Command: {cmd_str}

--- \uD83D\uDCCB Log Tail ---
{log_tail}
"""

print(f'   -> Python: {sys.version.split()[0]}')

# 1. PushPlus
if pushplus_token:
    try:
        data = urllib.parse.urlencode({
            'token': pushplus_token,
            'title': title,
            'content': msg.replace('\n', '<br>'), # HTML 换行
            'template': 'html'
        }).encode('utf-8')
        req = urllib.request.Request('http://www.pushplus.plus/send', data=data)
        urllib.request.urlopen(req)
        print('   -> PushPlus: Sent OK')
    except Exception as e:
        print(f'   -> PushPlus Error: {e}')

# 2. Bark
if bark_url:
    try:
        url = bark_url.rstrip('/')
        
        headers = {'Content-Type': 'application/json; charset=utf-8'}
        payload = {
            "title": title,
            "body": msg,
            # "group": "BioServer",  # 暂时注释掉分组，防止通知被折叠或误静音
            "icon": "https://cdn-icons-png.flaticon.com/512/103/103085.png",
            "level": "active",
            "sound": "minuet"
        }
        
        data = json.dumps(payload).encode('utf-8')
        req = urllib.request.Request(url, data=data, headers=headers, method='POST')
        
        urllib.request.urlopen(req)
        print('   -> Bark: Sent OK')
    except Exception as e:
        print(f'   -> Bark Error: {e}')

# 3. DingTalk
if dingtalk_url:
    try:
        headers = {'Content-Type': 'application/json'}
        data = json.dumps({
            'msgtype': 'text',
            'text': {'content': f'[Server] {title}\n\n{msg}'}
        }).encode('utf-8')
        req = urllib.request.Request(dingtalk_url, data=data, headers=headers)
        urllib.request.urlopen(req)
        print('   -> DingTalk: Sent OK')
    except Exception as e:
        print(f'   -> DingTalk Error: {e}')
EOF

# 设置 PYTHONUTF8=1 强制 Python 使用 UTF-8 处理文件 IO
export PYTHONUTF8=1

# 执行 Python 脚本，直接传递参数
$PY_CMD "$NOTIFY_SCRIPT" "$EXIT_CODE" "$TASK_NAME" "$TIME_STR" "$CMD_STR" "$LOG_TMP"

# 清理临时文件
rm -f "$LOG_TMP" "$NOTIFY_SCRIPT"
exit $EXIT_CODE
