import argparse

def main():
    parser = argparse.ArgumentParser(description="Description...")
    parser.add_argument("--name", type=str, default="World", help="Name to greet")
    args = parser.parse_args()
    name = args.name
    print("Hello", name)

if __name__ == "__main__":
    main()