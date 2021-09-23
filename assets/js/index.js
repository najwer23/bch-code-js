class BCH {

    // wielomiany sa w postaci 1 + x + x^2 + x^3 ...
    constructor (args) {
        this.galoisBase = 2;
        this.galoisPower = args.primitivePolynomialPeroid.toString(2).length;
        this.primitivePolynomialPeroid = this.galoisBase ** this.galoisPower - 1;
        this.primitivePolynomial = this.getPrimitivePolynomial();
        this.msg = args.msg.length < args.msgLength ? args.msg.padStart(args.msgLength, '0') : args.msg; // wiadomosc
        this.howManyErrors = args.howManyErrors; //zdolnosc korekcyjna
        this.controlPart = this.primitivePolynomialPeroid - this.msg.length
    }

    getPrimitivePolynomial() {
        let primitivePolyTest = "1".padStart(this.galoisPower+1, "0")

        primitivePolyTest = "1101"
        let primitivePolyTestPart = primitivePolyTest.substring(0,primitivePolyTest.length-1)
        let notZeroIndexPrimitivePolyTest = Array.from(primitivePolyTestPart, (c, i) => c === '1' ? i : -1).filter(i => i !== -1)
        
        let primitivePolynomialPeroid = 
        console.log(notZeroIndexPrimitivePolyTest)
        let peroidInf = "1".padEnd(this.galoisPower, "0")
        return primitivePolyTest
    }
};


//The load event fires when a given resource has loaded.
window.onload = () => {
    let objBCH = {
        primitivePolynomialPeroid: 2**3-1, //calkowoty mozliwy wektor kodowy
        msg: "11", // kodowana wiadomosc
        msgLength: 4, // calkowita mozliwa wiadomosc
        howManyErrors: 1, // liczby mozliwych bledow do skorygowania 
    }

    let bch = new BCH(objBCH)
    console.log(bch)
}

