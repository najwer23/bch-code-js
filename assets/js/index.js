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
                while ((i + c) < n && A[i + c][i] == 0)
                    c++;  

                if ((i + c) == n)
                    break;

                for (j = i, k = 0; k <= n; k++) {
                    let temp = A[j][k];
                    A[j][k] = A[j+c][k];
                    A[j+c][k] = temp;
                }
            }
    
            for (j = 0; j < n; j++) { 
                if (i != j) {
                    let p = this.customMod(A[j][i] / A[i][i],this.galoisBase);
                    
                    for (k = 0; k <= n; k++) {
                        A[j][k] = this.customMod(A[j][k] - this.customMod(A[i][k] * p,this.galoisBase),this.galoisBase);  
                    }              
                                
                }
            }
        }

        return A.map(c => c[A.length])
    }

    getMinimalPolys() {
        let minimalPolys = [];

        for (let k=0; k<this.alfasRoots.length; k++) {
            let nPoly = k; 
            let degreePoly = this.alfasRoots[nPoly].length;
            let A = [];

            // macierz kwadratowa z dolaczonym wektorem rozwiazan
            for (let i=degreePoly-1; i>0; i--) {
                A.push((this.alfasRoots[nPoly][0]*i) % this.primitivePolynomialPeroid)
            }
            A.push(0)
            for (let i=0; i<this.galoisPower - degreePoly; i++) {
                A.push(this.primitivePolynomialPeroid)
            }
            A.push(((this.alfasRoots[nPoly][0])*degreePoly) % this.primitivePolynomialPeroid)

            A = A.map(x=>this.alfas[x].split("")) 
         
            A = this.makeTranspose(A)
            A = this.runGaussGF(A)


            let minimalPoly = [1];
            for (let i=0; i<degreePoly; i++) {
                minimalPoly.push(A[i])
            }
            minimalPoly = minimalPoly.reverse()
            for (let i=0; i<(this.galoisPower - degreePoly); i++) {
                minimalPoly.push(0)
            }

            minimalPolys.push(minimalPoly)
        }

        for(let i=0; i<minimalPolys.length; i++) {
            let temp = minimalPolys[i];
            for(let j=i+1; j<minimalPolys.length; j++) {
                if (temp.join("") == minimalPolys[j].join("")) {
                    minimalPolys.splice(j,1);
                    j--;
                    continue;
                }
            }
        }

        return minimalPolys.map(x=>x.join(""))
    }

    getRootsOfMinimalPoly() {
        let alfasRoots = [];

        for (let i=1; i<this.primitivePolynomialPeroid; i++) {
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

        alfas.push("0".padStart(this.galoisPower, "0"))

        return alfas
    }

    getPossibleMsgLength() {
        return this.primitivePolynomialPeroid - this.galoisPower*this.howManyErrors;
    }

    getPrimitivePolynomial() {
        let primitivePolyTest = "1"+"1".padStart(this.galoisPower-1, "0")
        for (let i=0; ;i++) {
            if (this.checkIfPolynomialIsPrimitive(primitivePolyTest)) {
                console.log("dupa")
                break;
            }

            primitivePolyTest = primitivePolyTest.split("").map(x=>+x);
            primitivePolyTest[0] += 1; 

            for (let j=0; j<primitivePolyTest.length; j++) {
                if (primitivePolyTest[j]>=this.galoisBase) {
                    let r = this.customMod(primitivePolyTest[j], this.galoisBase)
                    primitivePolyTest[j+1] += 1
                    primitivePolyTest[j] = r
                }
            }

            primitivePolyTest = [...primitivePolyTest].join("")
        }

        return primitivePolyTest
    }

    checkIfPolynomialIsPrimitive(primitivePolyTest) {   
        // primitivePolyTest = "1000100001" 
        // primitivePolyTest = "100011101".split("").reverse().join("") 
        
        primitivePolyTest = primitivePolyTest.slice(0,this.galoisPower-1).split("")
        let indexs = primitivePolyTest.map((x,i)=>x==1?i:-1).filter(x=>x!==-1)

        let arr = "1".padEnd(this.galoisPower,"0").split("");
        for (let i=0; i<this.primitivePolynomialPeroid*2; i++) {
            let suma = 0;
            for (let j=0; j<indexs.length; j++) {
                suma = this.customMod(+arr[indexs[j]+i]+suma,this.galoisBase);
            }
            arr.push(suma)
        }

        this.cycleOfField = arr.slice(0,this.primitivePolynomialPeroid).join("")
        let roots = this.getElementsOfField()

        for (let i=0; i<roots.length; i++) {
            if ((roots[2] == roots[i]) && (i != 2))
                return false
        }
                
        return true
    }
};


//The load event fires when a given resource has loaded.
window.onload = () => {   
    let objBCH = {
        primitivePolynomialPeroid: 2**6-1, //calkowoty mozliwy wektor kodowy
        msg: "11", // kodowana wiadomosc
        howManyErrors: 1, // liczby mozliwych bledow do skorygowania 
    }

    let bch = new BCH(objBCH)
    console.log(bch)
    console.log(bch.getMinimalPolys())

}

