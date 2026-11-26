"""
LaTeX Auto-Builder (文件保存自动编译工具)
适用于 Antigravity IDE / VS Code / 终端

用法：
    python auto_build.py
当检测到 main.tex 或 reference.bib 保存修改时，会自动调用 xelatex 进行秒级增量编译。
"""

import os
import sys
import time
import subprocess

if sys.platform == 'win32':
    try:
        sys.stdout.reconfigure(encoding='utf-8')
    except Exception:
        pass

TEX_DIR = os.path.dirname(os.path.abspath(__file__))
MAIN_TEX = os.path.join(TEX_DIR, 'main.tex')
MAIN_PDF = os.path.join(TEX_DIR, 'main.pdf')
XELATEX_EXE = r'E:\texlive\2026\bin\windows\xelatex.exe'

def get_mtime(path):
    try:
        return os.path.getmtime(path)
    except OSError:
        return 0

def compile_tex():
    print(f"\n[{time.strftime('%H:%M:%S')}] 检测到文件保存，正在编译 main.pdf ...")
    start = time.time()
    cmd = [
        XELATEX_EXE,
        '-synctex=1',
        '-interaction=nonstopmode',
        '-file-line-error',
        'main.tex'
    ]
    try:
        proc = subprocess.run(cmd, cwd=TEX_DIR, capture_output=True, text=True, encoding='latin-1')
        duration = time.time() - start
        if proc.returncode == 0:
            print(f"[{time.strftime('%H:%M:%S')}] [OK] 编译成功！耗时 {duration:.2f} 秒 -> {MAIN_PDF}")
        else:
            print(f"[{time.strftime('%H:%M:%S')}] [ERROR] 编译报错 (退出码 {proc.returncode})，错误摘要：")
            for line in proc.stdout.splitlines()[-15:]:
                print("  ", line)
    except Exception as e:
        print(f"[{time.strftime('%H:%M:%S')}] 执行异常: {e}")

def main():
    print("=" * 60)
    print("  BioNLP LaTeX Auto-Builder 已启动")
    print(f"  监控目录: {TEX_DIR}")
    print(f"  监控主文件: main.tex")
    print(f"  编译器路径: {XELATEX_EXE}")
    print("  现在只需在 IDE 中按 Ctrl+S，这里就会自动重新编译！")
    print("  (按 Ctrl+C 可退出)")
    print("=" * 60)

    last_tex_mtime = get_mtime(MAIN_TEX)

    while True:
        try:
            time.sleep(0.5)
            curr_tex_mtime = get_mtime(MAIN_TEX)
            if curr_tex_mtime != last_tex_mtime and curr_tex_mtime > 0:
                last_tex_mtime = curr_tex_mtime
                # 稍微缓冲 0.2 秒以防编辑器正在写入
                time.sleep(0.2)
                compile_tex()
        except KeyboardInterrupt:
            print("\n监控已退出。")
            break

if __name__ == '__main__':
    main()
