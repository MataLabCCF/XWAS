def openListFile(fileName):
    listFile = open(fileName)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='PCA and regression')

    required = parser.add_argument_group("Required arguments")
    required.add_argument('-l', '--list', help='List of hybrid files. Format: <NAME> <HYBRID> <PGEN prefix file>',
                          required=True)
    required.add_argument('-f', '--folder', help='Output folder name', required=True)
    required.add_argument('-n', '--name', help='Name to output file', required=True)
    required.add_argument('-g', '--gwama', help='Path to gwama', required=True)
    required.add_argument('-p', '--plink2', help='Path to PLINK2', required=True)

    args = parser.parse_args()
    listFile = open(args.list)

    datasets = openListFile(args.list)