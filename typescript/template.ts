declare var require: (x: string) => any;

function Main(input: string[]) {


  console.log();
}

Main(require('fs').readFileSync('/dev/stdin', 'utf8').split('\n'));
