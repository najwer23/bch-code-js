class BCH {

    // wielomiany sa w postaci 1 + x + x^2 + x^3 ...
    constructor (args) {
        this.galoisBase = 2;
        this.galoisPower = args.primitivePolynomialPeroid.toString(2).length;
        this.primitivePolynomialPeroid = this.galoisBase ** this.galoisPower - 1;
        this.howManyErrors = args.howManyErrors; //zdolnosc korekcyjna
        this.msgLength = this.getPossibleMsgLength(); // calkowita mozliwa wiadomosc
        this.cycleOfField = "";
        this.primitivePolynomial = this.getPrimitivePolynomial();
        this.msg = args.msg.padStart(this.msgLength, '0'); // wiadomosc
        this.controlPart = this.primitivePolynomialPeroid - this.msg.length
        this.alfas = this.getElementsOfField()
        this.alfasRoots = this.getRootsOfMinimalPoly()
    }

    makeTranspose(arr) {
        return arr[0].map((x,i) => arr.map(x => x[i]));
    }

    getRootsOfMinimalPoly() {
        let alfasRoots = [];

        for (let i=1; i<this.primitivePolynomialPeroid/2; i=i+2) {
            let j=0;
            let alfaRootsOne = [];
            while(j<this.galoisPower) {
                alfaRootsOne.push((i*(this.galoisBase**(j))) % (this.primitivePolynomialPeroid))
                j++;
            }
            alfasRoots.push([...new Set(alfaRootsOne)])
        }

        return alfasRoots
    }

    getElementsOfField() {
        let alfas = [];
        let cycle = this.cycleOfField + this.cycleOfField

        for (let i=0; i<this.primitivePolynomialPeroid; i++) {
            alfas.push(cycle.slice(i,i+this.galoisPower))
        }

        return alfas
    }

    getPossibleMsgLength() {
        return this.primitivePolynomialPeroid - this.galoisPower*this.howManyErrors;
    }

    getPrimitivePolynomial() {
        let primitivePolyTest = "1".padStart(this.galoisPower+1, "0")
        for (let i=0; ;i++) {
            if (this.checkIfPolynomialIsPrimitive(primitivePolyTest)) {
                break;
            }

            primitivePolyTest = primitivePolyTest.split("").map(x=>+x);
            primitivePolyTest[0] += 1; 

            for (let j=0; j<primitivePolyTest.length; j++) {
                if (primitivePolyTest[j]>=this.galoisBase) {
                    let r = primitivePolyTest[j] % this.galoisBase
                    primitivePolyTest[j+1] += 1
                    primitivePolyTest[j] = r
                }
            }

            primitivePolyTest = [...primitivePolyTest].join("")
        }

        return primitivePolyTest
    }

    checkIfPolynomialIsPrimitive(primitivePolyTest) {
        let primitivePolyTestPart = primitivePolyTest.substring(0,primitivePolyTest.length-1)
        let notZeroIndexPrimitivePolyTest = Array.from(primitivePolyTestPart, (c, i) => c === '1' ? i : -1).filter(i => i !== -1)

        let arrCyclicPeroid = "1".padEnd(this.galoisPower, "0").split("")
        for (let i=0; i<this.primitivePolynomialPeroid*2; i++) {
            let sum = notZeroIndexPrimitivePolyTest.reduce((a,b)=>(a + +arrCyclicPeroid[b+i]) % this.galoisBase,0)
            arrCyclicPeroid.push(sum)
        }

        this.cycleOfField = arrCyclicPeroid.slice(0,this.primitivePolynomialPeroid).join("")
        return (arrCyclicPeroid.splice(0,this.primitivePolynomialPeroid).join("") == arrCyclicPeroid.splice(0,this.primitivePolynomialPeroid).join(""))
    }
};


//The load event fires when a given resource has loaded.
window.onload = () => {
    let objBCH = {
        primitivePolynomialPeroid: 2**4-1, //calkowoty mozliwy wektor kodowy
        msg: "11", // kodowana wiadomosc
        howManyErrors: 1, // liczby mozliwych bledow do skorygowania 
    }

    let bch = new BCH(objBCH)
    console.log(bch)
}

