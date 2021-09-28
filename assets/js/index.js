const cloneDeep = require('lodash.clonedeep')

class BCH {

    // wielomiany sa w postaci 1 + x + x^2 + x^3 ...
    constructor (args) {

        let date1 = new Date();

        this.galoisBase = 2;
        this.galoisPower = args.codeLength.toString(2).length;
        this.codeLength = this.galoisBase ** this.galoisPower - 1;
        this.howManyErrors = args.howManyErrors; //zdolnosc korekcyjna
        this.cycleFromPseudoRandomGenerator = "";

        let date7 = new Date();
        this.primitivePolynomial = this.getPrimitivePolynomial();
        let date8 = new Date();
        this.timePrimitivePolynomial = (date8 - date7)/1000

        this.alfas = this.getElementsOfField()
        this.minimalPolynomialsRootsAsAlfa = this.getRootsOfMinimalPoly()

        let date9 = new Date();
        this.minimalPolynomials = this.getMinimalPolynomials()
        let date10 = new Date();
        this.timeMinimalPolynomials = (date10 - date9)/1000

        this.polynomialGeneratingCode = this.getPolynomialGeneratingCode()
        this.controlPart = this.polynomialGeneratingCode.length
        this.msgLength = this.codeLength - this.controlPart + 2// calkowita mozliwa wiadomosc - stopien wielominau, jesli wielomian jest 8 stopnia to jest 9 bitow
        this.msgPaddingAtStart = this.msgLength - args.msg.length
        this.msg = args.msg.padStart(this.msgLength, '0'); // wiadomosc
        this.msgWithoutPadding = args.msg; // wiadomosc

        let date3 = new Date();
        this.encodeMSG = this.encodeMsgBCH()
        let date4 = new Date();
        this.timeEncodeMsgBCH = (date4 - date3)/1000
        
        let date5 = new Date();
        this.decodeMSG = this.decodeMsgBCH()
        let date6 = new Date();
        this.timeDecodeMsgBCH = (date6 - date5)/1000

        this.syndrom;

        let date2 = new Date();
        this.timeAllBCH = (date2 - date1)/1000;

        this.timeMulDecodeMsgBCH;
        this.timeDivDecodeMsgBCH;

        this.showTime();
    }

    showTime() {
        console.log("BCH - short")
        console.log("Czas kodowanie: ",  this.timeEncodeMsgBCH)
        console.log("Czas znalezienia wiel. prym: ",  this.timePrimitivePolynomial)
        console.log("Czas znalezienia wiel. minimalnych: ",  this.timeMinimalPolynomials)
        console.log("Czas mnozenie msg*x^m:",  this.timeMulDecodeMsgBCH)
        console.log("Czas dzielienia msg*x^m / G(x): ",  this.timeDivDecodeMsgBCH)
        console.log("Czas dekodowania: ",  this.timeDecodeMsgBCH)
        console.log("Czas całkowity: ",  this.timeAllBCH)
    }

    encodeMsgBCH() {
        //wielomian stopnia generujacego 
        let X = "1".padStart(this.controlPart, "0")

        let date1 = new Date();
        let Y = this.mul2Polynomials(this.msg, X)
        let date2 = new Date();
        this.timeMulDecodeMsgBCH = (date2 - date1)/1000

        let date3 = new Date()
        let Z = this.div2Polynomials(Y,this.polynomialGeneratingCode)
        let date4 = new Date()
        this.timeDivDecodeMsgBCH = (date4 - date3)/1000

        return this.add2Polynomials(Y,Z.remainder)
    }

    createMatrixOfSyndroms() {
        let s = this.syndrom.split("");
        let arrS = [];
        let arrSRow = []
        let arrSsingleSyndrom = []

        for (let i=0; i<this.howManyErrors; i++) {
            arrSRow = [];
            for(let j=0; j<this.howManyErrors + 1; j++) {
                arrSsingleSyndrom = []
                //przeksztalc syndrom na alfe, ale na wielomian
                for (let k=0; k<s.length; k++) {
                    if (s[k]=="1") {
                        let S = (+(k)*(j+1+i)) % this.codeLength
                        S = this.alfas[S]
                        arrSsingleSyndrom.push(S)
                    }
                }

                // superpozycja elemetow ciala - alfy i wilomian
                arrSsingleSyndrom = arrSsingleSyndrom.reduce((a,b)=>a = this.add2Polynomials(a,b))

                arrSsingleSyndrom = {
                    poly: arrSsingleSyndrom,
                    alfa: this.alfas.indexOf(arrSsingleSyndrom),
                }
                
                arrSRow.push(arrSsingleSyndrom)
            }
            arrS.push(arrSRow)
        }

        return arrS
    }

    //gauss jordan elimination dla wyznacznika i rozwiazania ukladu rownan w ciele pod funkcje lambda
    runGaussGFForPolynomial(matrixOfSyndorms) {
        let resultOfDivAlfas;
        let resultOfMulAlfas;
        let c=0;
        let A;
        let temp;

        for (let t = this.howManyErrors+1; t>1; t--) { 
            A = matrixOfSyndorms.slice(0, t-1).map(i => i.slice(0, t));
            c=0;
            for (let i = 0; i < A.length; i++) {
                if (A[i][i].alfa == this.codeLength) {
                    c = 1;
                    
                    while ((i + c) < A.length && A[i + c][i].alfa == this.codeLength) { c++; }
                    
                    if ((i + c) == A.length) { break; }
                    
                    //pivot
                    for (j = i, k = 0; k <= A.length; k++) {
                        temp = A[j][k];
                        A[j][k] = A[j+c][k];
                        A[j+c][k] = temp;
                    }
                }
        
                for (let j = 0; j < A.length; j++) { 
                    if (i != j) {
                        resultOfDivAlfas = this.operationOn2Alfas(A[j][i], A[i][i], "/") 
                        for (let k = 0; k <= A.length; k++) {
                            resultOfMulAlfas = this.operationOn2Alfas(resultOfDivAlfas,A[i][k], "*") 
                            A[j][k] = this.operationOn2Alfas(A[j][k], resultOfMulAlfas, "+")
                        }
                    }
                }
            }

            //determinant(A)
            let det = "1";
            for (let i=0; i<t-1; i++) {
                det = this.mul2Polynomials(det,A[i][i].poly)
            }

            if(det.lastIndexOf("1")>0) {
                //lambdaCoefficients
                for (let i=0; i<A.length; i++) {
                    A[i][A.length] = this.operationOn2Alfas(A[i][A.length], A[i][i],"/")
                }
                return A.map(c => c[A.length])
            }
        }
    }

    makeOppositeBit(s, index) {
        return s.substring(0, index) + (s[index] == "1" ? "0" : "1") + s.substring(index + 1);
    }

    decodeMsgBCH() {
        //testc

        //make errors in msg
        this.encodeMSG = this.makeOppositeBit(this.encodeMSG,1)
        this.encodeMSG = this.makeOppositeBit(this.encodeMSG,2)
        this.encodeMSG = this.makeOppositeBit(this.encodeMSG,3)
        
        let Cy = this.div2Polynomials(this.encodeMSG,this.polynomialGeneratingCode)
        this.syndrom = Cy.remainder

        //brak bledu w wektorze kodowym
        if (this.syndrom == 0) {
            return this.encodeMSG.slice(this.encodeMSG.length-this.msgLength+this.msgPaddingAtStart,this.encodeMSG.length)
        }

        let matrixOfSyndorms = this.createMatrixOfSyndroms(this.syndrom)
        //matrixOfSyndorms.map(x=>console.log("Macierz syndromow", x.map(y=>y.alfa)))
      


        //console.log("Lambda")
        let lambdaCoefficients = this.runGaussGFForPolynomial(matrixOfSyndorms)
        lambdaCoefficients.push({ poly: this.alfas[0], alfa: 0})
        lambdaCoefficients.reverse();

        //1+((12+8)%15)*x+((9+16)%15)*x^2
        for (let i=0; i<=this.codeLength-1; i++) {
            let p = "0"
            for (let j=0; j<lambdaCoefficients.length; j++) {
                let s = (lambdaCoefficients[j].alfa+(i*(j))) % this.codeLength
                p = this.add2Polynomials(p,this.alfas[s])
                //console.log(i,j,s, this.alfas[s], lambdaCoefficients[j].alfa)
            }
            if (p.lastIndexOf("1") < 0) {
                console.log("Blad w pozycji", this.customMod(-1*i,this.codeLength))
            }
        }
       
        
    }

    operationOn2Alfas(a,b,operation) {
        let n = this.codeLength;
        let r_alfa;

        if (operation == "+") {
            let resultPoly = this.add2Polynomials(a.poly,b.poly)
            return {
                poly: resultPoly,
                alfa: this.alfas.indexOf(resultPoly)
            }
        }

        if (operation == "/") {
            r_alfa = (a.alfa == n || b.alfa == n) ? n : this.customMod((a.alfa - b.alfa), n) 
        }

        if (operation == "*") {
            r_alfa = (a.alfa == n || b.alfa == n) ? n : this.customMod((a.alfa + b.alfa), n) 
        }

        return {
            poly: this.alfas[r_alfa],
            alfa: r_alfa
        }
    }


    getPolynomialGeneratingCode() {
        let G = this.minimalPolynomials[0]
        for (let i=1; i<this.howManyErrors; i++) {
            G = this.mul2Polynomials(G,this.minimalPolynomials[i])
        }
        return G.slice(0,G.lastIndexOf("1")+1)
    }

    makeTranspose(arr) {
        return arr[0].map((x,i) => arr.map(x => x[i]));
    }

    // JavaScript ma problem z operacja modula na liczbach ujemnych
    customMod(a,n) {
        return ((a%n)+n)%n;
    }

    // Gauss-Jordan Elimination Method in GF
    runGaussGF(A) {
        let i=0;
        let j=0;
        let k=0;
        let c=0;
        let n = A.length;

        for (i = 0; i < n; i++) {
            if (A[i][i] == 0) {
                c = 1;
                while ((i + c) < n && A[i + c][i] == 0) c++;  
                
                if ((i + c) == n) break;

                for (j = i, k = 0; k <= n; k++) {
                    let temp = A[j][k];
                    A[j][k] = A[j+c][k];
                    A[j+c][k] = temp;
                }
            }
    
            for (j = 0; j < n; j++) { 
                if (i != j) {
                    let p = this.customMod(A[j][i] / A[i][i],this.galoisBase);
                    for (k = 0; k <= n; k++) 
                        A[j][k] = this.customMod(A[j][k] - this.customMod(A[i][k] * p,this.galoisBase),this.galoisBase);                     
                }
            }
        }

        return A.map(c => c[A.length])
    }

    getMinimalPolynomials() {
        let minPolynomials = [];

        for (let k=0; k<this.minimalPolynomialsRootsAsAlfa.length; k++) {
            let degreePolynomial = this.minimalPolynomialsRootsAsAlfa[k].length;
            let A = [];

            // macierz kwadratowa wektorow alf
            for (let i=degreePolynomial-1; i>0; i--) {
                A.push((this.minimalPolynomialsRootsAsAlfa[k][0]*i) % this.codeLength)
            }
            A.push(0)
            for (let i=0; i<this.galoisPower - degreePolynomial; i++) {
                A.push(this.codeLength)
            }

            // wektor rozwiazan 
            A.push(((this.minimalPolynomialsRootsAsAlfa[k][0])*degreePolynomial) % this.codeLength)

            //znajdz wspolczyniki wielomianu minimanego - rozwiazanie liniowego ukladu rownan 
            A = A.map(x=>this.alfas[x].split("")) 
            A = this.makeTranspose(A)
            A = this.runGaussGF(A)

            // transformuj wektor rozwiazan do postaci 1 + x + x^2...
            let minPolynomial = [1];
            for (let i=0; i<degreePolynomial; i++) {
                minPolynomial.push(A[i])
            }
            minPolynomial = minPolynomial.reverse()
            for (let i=0; i<(this.galoisPower - degreePolynomial); i++) {
                minPolynomial.push(0)
            }
            
            //dodaj wilomian minimalny do tablicy
            minPolynomials.push(minPolynomial)
        }

        //usun duplikaty 
        for(let i=0; i<minPolynomials.length; i++) {
            let temp = minPolynomials[i];
            for(let j=i+1; j<minPolynomials.length; j++) {
                if (temp.join("") == minPolynomials[j].join("")) {
                    minPolynomials.splice(j,1); j--;
                    continue
                }
            }
        }

        //1 + x + x^2..
        return minPolynomials.map(x=>x.join(""))
    }

    getRootsOfMinimalPoly() {
        let minimalPolynomialsRootsAsAlfa = [];

        for (let i=1; i<this.codeLength; i++) {
            let j=0;
            let alfaRootsOne = [];
            while(j<this.galoisPower) {
                alfaRootsOne.push((i*(this.galoisBase**(j))) % (this.codeLength))
                j++;
            }
            minimalPolynomialsRootsAsAlfa.push([...new Set(alfaRootsOne)])
        }

        return minimalPolynomialsRootsAsAlfa
    }

    getElementsOfField() {
        let alfas = [];
        let cycle = this.cycleFromPseudoRandomGenerator + this.cycleFromPseudoRandomGenerator

        for (let i=0; i<this.codeLength; i++) {
            alfas.push(cycle.slice(i,i+this.galoisPower))
        }

        alfas.push("0".padStart(this.galoisPower, "0"))

        return alfas
    }

    getPrimitivePolynomial() {
        //poczatkowy wielomian 
        let primitivePolynomialTest = "1"+"1".padStart(this.galoisPower, "0")

        //szukaj wielomianu prymitywnego tak dlugo az znajdziesz
        for (let i=0; ;i++) {
            if (this.checkIfPolynomialIsPrimitive(primitivePolynomialTest)) {
                break;
            }

            primitivePolynomialTest = primitivePolynomialTest.split("").map(x=>+x);
            primitivePolynomialTest[0] += 1; 

            // inkrementowanie po stopniach wielomianu, jednak dalej w obrebie ciala
            for (let j=0; j<primitivePolynomialTest.length; j++) {
                if (primitivePolynomialTest[j]>=this.galoisBase) {
                    let r = this.customMod(primitivePolynomialTest[j], this.galoisBase)
                    primitivePolynomialTest[j+1] += 1
                    primitivePolynomialTest[j] = r
                }
            }
            primitivePolynomialTest = [...primitivePolynomialTest].join("")
        }
        return primitivePolynomialTest
    }

    checkIfPolynomialIsPrimitive(primitivePolynomialTest) {   
        //wytnij najwyzsza potege wielominau
        primitivePolynomialTest = primitivePolynomialTest.slice(0,this.galoisPower-1).split("")
        
        //pseudo generator losowy elementow ciala
        let indexs = primitivePolynomialTest.map((x,i)=>x==1?i:-1).filter(x=>x!==-1)
        let arr = "1".padEnd(this.galoisPower,"0").split("");
        for (let i=0; i<this.codeLength; i++) {
            let suma = 0;
            for (let j=0; j<indexs.length; j++) {
                suma = this.customMod(+arr[indexs[j]+i]+suma,this.galoisBase);
            }
            arr.push(suma)
        }
        this.cycleFromPseudoRandomGenerator = arr.slice(0,this.codeLength).join("")
       
        //sprawdz czy utworzone elementy z losowego generatora sa unikatowe
        //sprawdz czy wielomian prymitywny utworzyl maksymalny okres
        let roots = this.getElementsOfField()
        for (let i=0; i<roots.length; i++) {
            if ((roots[2] == roots[i]) && (i != 2)) {
                return false
            }
        }      
        return true
    }

    shiftStringRight(s, shift) {
        let l = s.length - shift;
        return s.substring(l) + s.substring(0, l)
    }

    getPolynomialDegreeDifference(str1, str2) {
        return str1.lastIndexOf('1') - str2.lastIndexOf('1')
    }

    div2Polynomials(a,b) {
        b = b.padEnd(a.length,"0")
        let copyb = b;
        let isDive = true;
        let polynomialDegreeDifference = this.getPolynomialDegreeDifference(a,b)
        let r = Array(Math.abs(polynomialDegreeDifference+1)).fill(0);

        while (isDive) {
            b = copyb;
            polynomialDegreeDifference = this.getPolynomialDegreeDifference(a,b)
            if (polynomialDegreeDifference >= 0) {
                b = this.shiftStringRight(b,polynomialDegreeDifference) 
                r[polynomialDegreeDifference]=1
                a = this.add2Polynomials(a,b)
            }
            isDive = !(polynomialDegreeDifference < 0)  
        }

        r = r.join("")

        return {
            "result": r.slice(0,r.lastIndexOf("1")+1),
            "remainder": a.slice(0,a.lastIndexOf("1")+1) || "0"
        }
    }

    mul2Polynomials(a,b) {
        a = a.split("")
        b = b.split("")
        let r;
        let R = [];

        let k=0;
        for (let i=0; i<a.length; i++) {
            if (a[i] != "0") {
                r = "0".padEnd(a.length+b.length-1, '0').split("");
                for (let j=0; j<b.length; j++) {
                    r[j+k]=a[i]*b[j];
                }
                R.push(r.join(""))
            }
            k++;   
        }  

        let result = "0"
        for (let i=0; i<R.length; i++) {
            result = this.add2Polynomials(result,R[i])
        }
        
        return result
    }

    add2Polynomials(a,b) {
        let r = a.length >= b.length ? a : b;
        a = a.padEnd(r.length,'0')
        b = b.padEnd(r.length,'0')
        r = r.split("");
        for (let i=0; i<r.length; i++) r[i] = a[i] ^ b[i];
        return r.join("");
    }
};



window.onload = () => {   
    let objBCH = {
        codeLength: 2**14-1, //calkowoty mozliwy wektor kodowy
        msg: "111", // kodowana wiadomosc
        howManyErrors: 3, // liczby mozliwych bledow do skorygowania 
    }

    let bch = new BCH(objBCH)
    console.log(bch)
}

