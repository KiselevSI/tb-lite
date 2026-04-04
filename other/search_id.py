#!/usr/bin/env python3
import argparse
import os
import re
from pathlib import Path


def read_ids(file_path: str) -> list[str]:
    ids = []
    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            value = line.strip()
            if value:
                ids.append(value)
    return ids


def build_pattern(ids: list[str]) -> re.Pattern:
    escaped = [re.escape(x) for x in ids]
    return re.compile("|".join(escaped))


def main():
    parser = argparse.ArgumentParser(
        description="Поиск файлов и папок по списку ID в именах"
    )
    parser.add_argument("-i", required=True, help="Входная папка")
    parser.add_argument("-l", required=True, help="TXT файл со списком ID")
    parser.add_argument("-o", required=True, help="Выходной TXT файл")
    args = parser.parse_args()

    input_dir = Path(args.i)
    ids_file = Path(args.l)
    output_file = Path(args.o)

    if not input_dir.is_dir():
        print(f"Ошибка: входная папка не существует или это не папка: {input_dir}")
        return

    if not ids_file.is_file():
        print(f"Ошибка: файл со списком ID не существует: {ids_file}")
        return

    ids = read_ids(str(ids_file))
    ids = list(dict.fromkeys(ids))

    if not ids:
        print("Ошибка: список ID пуст")
        return

    pattern = build_pattern(ids)

    found_paths = []
    found_set = set()

    for root, dirs, files in os.walk(input_dir, topdown=True):
        matched_dirs = []

        for dirname in dirs:
            if pattern.search(dirname):
                dirpath = os.path.normpath(os.path.join(root, dirname))
                if dirpath not in found_set:
                    found_paths.append(dirpath)
                    found_set.add(dirpath)
                matched_dirs.append(dirname)

        dirs[:] = [d for d in dirs if d not in matched_dirs]

        for filename in files:
            if pattern.search(filename):
                filepath = os.path.normpath(os.path.join(root, filename))
                if filepath not in found_set:
                    found_paths.append(filepath)
                    found_set.add(filepath)

    with open(output_file, "w", encoding="utf-8") as f:
        for path in found_paths:
            f.write(path + "\n")

    print(f"Готово. Найдено путей: {len(found_paths)}")
    print(f"Результат записан в: {output_file}")


if __name__ == "__main__":
    main()