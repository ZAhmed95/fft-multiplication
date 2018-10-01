import java.awt.*;
import java.awt.event.*;
import java.math.BigInteger;
import javax.swing.*;
import javax.swing.table.*;

public class FFTMultiplicationMod {
	
	public static void main(String[] args)
	{
		SwingUtilities.invokeLater(new Runnable()
		{
			@Override
			public void run()
			{
				initGUI();
			}
		});
	}
	
	/*
	 * =========================================================================================
	 * FFT CODE
	 * =========================================================================================
	 */
	//FFT multiplication algorithm to multiply a and b
	static BigInteger multiply(BigInteger a, BigInteger b)
	{
		//get signs of both inputs
		int signA = a.signum();
		int signB = b.signum();
		//calculate sign of product
		int signC = signA * signB;
		//convert inputs to their absolute values
		a = a.multiply(BigInteger.valueOf(signA));
		b = b.multiply(BigInteger.valueOf(signB));
		
		int digitsA = digitCount(a, BigInteger.TEN);
		int digitsB = digitCount(b, BigInteger.TEN);
		//trivial cases:
		if (digitsA <= 3 && digitsB <= 3) 
		{
			fftTable = null;
			return a.multiply(b);
		}
		
		//K must be large enough to hold multiplication of a and b without losing any information
		int K = nextPowOf2(digitsA + digitsB);
		//mod = 2^K + 1
		BigInteger mod = BigInteger.ONE.shiftLeft(K).add(BigInteger.ONE);
		
		//convert integers a and b into convolution vectors
		BigInteger[] vectorA = formatNumber(a, K);
		BigInteger[] vectorB = formatNumber(b, K);
		
		//because the modulo, call it p, is equal to 2^K + 1, 2^k = -1 (mod p), which means 2^(2K) = 4^K = 1 (mod p)
		//so our primitive Kth root of unity mod p is 4.
		BigInteger w = BigInteger.valueOf(4);
		//calculate the FFTs of a and b
		BigInteger[] f_A = FFT(vectorA, K, w, mod);
		BigInteger[] f_B = FFT(vectorB, K, w, mod);
		
		//get the convolution f_C = f_A * f_B by pointwise multiplication of elements of f_A and f_B
		BigInteger[] f_C = new BigInteger[K];
		for (int i = 0; i < K; i++) f_C[i] = f_A[i].multiply(f_B[i]).mod(mod);
		
		//do the inverse FFT on f_C to get vectorC
		//inverse FFT is in fact 1/K * FFT on f_C using w^-1
		//Inverse of K is found using java's built in BigInteger.modInverse() method
		BigInteger K_inv = BigInteger.valueOf(K).modInverse(mod);
		
		//inverse of 4 in mod p is 4^(K-1), or 2^(2K - 2), which is why we put in 1 << (2*K - 2)
		BigInteger[] vectorC = FFT(f_C, K, BigInteger.ONE.shiftLeft(2*K - 2), mod);
		
		for (int i = 0; i < K; i++)
		{
			vectorC[i] = vectorC[i].multiply(K_inv).mod(mod);
		}
		
		//obtain product "c" by evaluating the polynomial defined by vectorC at x = 10
		BigInteger c = BigInteger.ZERO;
		BigInteger x = BigInteger.ONE;
		for (int i = 0; i < K; i++)
		{
			c = c.add(x.multiply(vectorC[i]));
			x = x.multiply(BigInteger.TEN);
		}
		
		//create 2D array of the data to post to fftTable
		BigInteger[][] data = new BigInteger[K][6];
		for (int i = 0; i < K; i++)
		{
			data[i][0] = vectorA[i]; //vector representation of digits of A
			data[i][1] = vectorB[i]; //vector representation of digits of B
			data[i][2] = f_A[i]; //FFT of A
			data[i][3] = f_B[i]; //FFT of B
			data[i][4] = f_C[i]; //FFT of C (calculated by F(C) = F(A) * F(B))
			data[i][5] = vectorC[i]; //calculated by IFFT(F(C))
		}
		updateTable(data);
		
		//return c multiplied by its sign
		return c.multiply(BigInteger.valueOf(signC));
	}
	
	//FFT algorithm, calculates the DFT of input vector a
	static BigInteger[] FFT(BigInteger[] a, int K, BigInteger w, BigInteger mod)
	{
		//base case
		if (K == 1)
		{
			return a;
		}
		//Construct the even and odd vectors of a:
		BigInteger[] aEven = new BigInteger[K/2];
		BigInteger[] aOdd = new BigInteger[K/2];
		for (int i = 0; i < K/2; i++)
		{
			aEven[i] = a[2 * i];
			aOdd[i] = a[2 * i + 1];
		}
		//recursively evaluate the FFT of each half, even and odd.
		BigInteger[] aEvenFFT = FFT(aEven, K/2, w.multiply(w).mod(mod), mod);
		BigInteger[] aOddFFT = FFT(aOdd, K/2, w.multiply(w).mod(mod), mod);
		BigInteger[] aFFT = new BigInteger[K];
		//recombine them to get our desired DFT, using the equation: A(x) = AEven(x^2) + x*AOdd(x^2)
		BigInteger x = BigInteger.ONE;
		for (int i = 0; i < K/2; i++)
		{
			//A(x) = AEven(x) + x * AOdd(x)
			aFFT[i] = (aEvenFFT[i].add(x.multiply(aOddFFT[i]))).mod(mod);
			//Due to the symmetry of the roots of unity, the second half of the evaluation points
			//are the negative of the first half, and are evaluated as A(-x) = AEven(x^2) - x * AOdd(x^2)
			aFFT[i + K/2] = (aEvenFFT[i].subtract(x.multiply(aOddFFT[i]))).mod(mod).add(mod).mod(mod);
			//update x
			x = x.multiply(w).mod(mod);
		}
		return aFFT;
	}
	/*
	 * =========================================================================================
	 * FFT CODE END
	 * =========================================================================================
	 */
	
	/*
	 * =========================================================================================
	 * HELPER METHODS
	 * =========================================================================================
	 */
	//this method calculates the next power of 2 that is greater than input integer a
	static int nextPowOf2(int a)
	{
		int b = 1;
		while (b < a) b = b << 1;
		return b;
	}
	
	//this method converts the given integer into a vector of digits with m elements.
	//m can be greater than the number of digits in b, the rest will just be padded zeroes
	static BigInteger[] formatNumber(BigInteger b, int m)
	{
		BigInteger[] vector = new BigInteger[m];
		for (int i = 0; i < m; i++)
		{
			vector[i] = b.mod(BigInteger.TEN);
			b = b.divide(BigInteger.TEN);
		}
		return vector;
	}
	
	//counts how many digits n has in base b
	static int digitCount(BigInteger n, BigInteger b)
	{
		int count = 0;
		while (n.compareTo(BigInteger.ZERO) > 0)
		{
			n = n.divide(b);
			count++;
		}
		return count;
	}
	
	//calculates log base 2 of K (rounded down)
	static int log2(int K)
	{
		int log = 0;
		while (K > 1)
		{
			K = K >> 1;
			log++;
		}
		return log;
	}
	/*
	 * =========================================================================================
	 * HELPER METHODS END
	 * =========================================================================================
	 */
	
	/*
	 * =========================================================================================
	 * GUI CODE
	 * =========================================================================================
	 */
	private static JFrame gui;
	private static JTextField txtA;
	private static JTextField txtB;
	private static JTextField txtC;
	private static JTable fftTable;
	//fftTable column names
	private static String[] columns = {"A", "B", "FFT(A)", "FFT(B)", "FFT(C)", "C"};
	private static JScrollPane jsp;
	
	static void initGUI()
	{
		gui = new JFrame();
		gui.setLayout(new FlowLayout());
		gui.setSize(750, 550);
		gui.setResizable(false);
		gui.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		//Label A
		JLabel lblA = new JLabel("A: ");
		gui.add(lblA);
		
		//the action listener used for user input of A and B
		MultiplyAction ma = new FFTMultiplicationMod().new MultiplyAction();
		
		//Textfield A
		txtA = new JTextField(15);
		txtA.addActionListener(ma);
		gui.add(txtA);
		
		//Label B
		JLabel lblB = new JLabel("B: ");
		gui.add(lblB);
		
		//Textfield B
		txtB = new JTextField(15);
		txtB.addActionListener(ma);
		gui.add(txtB);
		
		//Label A * B =
		JLabel lblAB = new JLabel("A * B = C: ");
		gui.add(lblAB);
		
		//Textfield for output of A * B
		txtC = new JTextField(15);
		gui.add(txtC);
		
		//Table to show data from FFT
		TableModel tbm = new DefaultTableModel(columns, 0);
		fftTable = new JTable(tbm);
		
		//scroll pane for table
		jsp = new JScrollPane();
		jsp.setViewportView(fftTable);
		gui.add(jsp);
				
		gui.setVisible(true);
	}
	
	//method called when multiplying A * B, updates the table to show vector form of A and B,
	//the FFT of A and B, the resulting FFT of C, and finally the vector form of C obtained after the inverse FFT
	static void updateTable(BigInteger[][] data)
	{
		fftTable = new JTable(data, columns);
		jsp.setViewportView(fftTable);
	}
	
	class MultiplyAction implements ActionListener
	{
		@Override
		public void actionPerformed(ActionEvent e)
		{
			BigInteger a = null;
			BigInteger b = null;
			BigInteger c = null;
			try
			{
				a = new BigInteger(txtA.getText());
				b = new BigInteger(txtB.getText());
			}
			catch(NumberFormatException nfe)
			{
				return;
			}
			
			c = multiply(a, b);
			
			txtC.setText(c.toString());
		}
	}
	/*
	 * =========================================================================================
	 * GUI CODE END
	 * =========================================================================================
	 */
}

