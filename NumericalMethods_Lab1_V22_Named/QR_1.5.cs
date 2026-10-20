using System.Numerics;

public static class QR
{
    public static (double[,] Q, double[,] R) Decompose(double[,] source)
    {
        Matrix.Square(source);
        int n = source.GetLength(0);
        double[,] q = Matrix.Identity(n), r = (double[,])source.Clone();
        for (int k = 0; k < n - 1; k++)
        {
            double norm = 0;
            for (int i = k; i < n; i++) norm = Math.Sqrt(norm * norm + r[i, k] * r[i, k]);
            if (norm < 1e-14) continue;
            double[] v = new double[n];
            v[k] = r[k, k] + (r[k, k] >= 0 ? norm : -norm);
            for (int i = k + 1; i < n; i++) v[i] = r[i, k];
            double length = 0;
            for (int i = k; i < n; i++) length += v[i] * v[i];
            length = Math.Sqrt(length);
            if (length < 1e-14) continue;
            for (int i = k; i < n; i++) v[i] /= length;
            for (int j = k; j < n; j++)
            {
                double dot = 0;
                for (int i = k; i < n; i++) dot += v[i] * r[i, j];
                for (int i = k; i < n; i++) r[i, j] -= 2 * v[i] * dot;
            }
            for (int i = 0; i < n; i++)
            {
                double dot = 0;
                for (int j = k; j < n; j++) dot += q[i, j] * v[j];
                for (int j = k; j < n; j++) q[i, j] -= 2 * dot * v[j];
            }
        }
        return (q, r);
    }

    public static (Complex[] values, int steps) Eigenvalues(double[,] source, double eps, int maxSteps)
    {
        Matrix.Square(source);
        int n = source.GetLength(0), m = n, steps = 0;
        double[,] a = (double[,])source.Clone();
        Complex[] values = new Complex[n];
        while (m > 0)
        {
            if (m == 1) { values[0] = a[0, 0]; break; }
            if (Small(a[m - 1, m - 2], a[m - 2, m - 2], a[m - 1, m - 1], eps))
            {
                values[m - 1] = a[m - 1, m - 1];
                m--;
                continue;
            }
            if (m == 2 || Small(a[m - 2, m - 3], a[m - 3, m - 3], a[m - 2, m - 2], eps))
            {
                (values[m - 2], values[m - 1]) = Pair(a[m - 2, m - 2], a[m - 2, m - 1], a[m - 1, m - 2], a[m - 1, m - 1]);
                m -= 2;
                continue;
            }
            if (steps >= maxSteps) throw new InvalidOperationException("QR-итерации не сошлись");
            double a11 = a[m - 2, m - 2], a12 = a[m - 2, m - 1];
            double a21 = a[m - 1, m - 2], a22 = a[m - 1, m - 1];
            double trace = a11 + a22, determinant = a11 * a22 - a12 * a21;
            double discriminant = trace * trace - 4 * determinant;
            double[,] block = new double[m, m];
            for (int i = 0; i < m; i++)
                for (int j = 0; j < m; j++) block[i, j] = a[i, j];
            double[,] next;
            if (discriminant >= 0)
            {
                double root = Math.Sqrt(discriminant);
                double l1 = (trace + root) / 2, l2 = (trace - root) / 2;
                double shift = Math.Abs(l1 - a22) < Math.Abs(l2 - a22) ? l1 : l2;
                for (int i = 0; i < m; i++) block[i, i] -= shift;
                var (q, r) = Decompose(block);
                next = Matrix.Multiply(r, q);
                for (int i = 0; i < m; i++) next[i, i] += shift;
            }
            else
            {
                double[,] polynomial = Matrix.Multiply(block, block);
                for (int i = 0; i < m; i++)
                {
                    for (int j = 0; j < m; j++) polynomial[i, j] -= trace * block[i, j];
                    polynomial[i, i] += determinant;
                }
                var (q, _) = Decompose(polynomial);
                next = Matrix.Multiply(Matrix.Multiply(Matrix.Transpose(q), block), q);
            }
            for (int i = 0; i < m; i++)
                for (int j = 0; j < m; j++) a[i, j] = next[i, j];
            steps++;
        }
        return (values, steps);
    }

    static bool Small(double value, double a, double b, double eps)
        => Math.Abs(value) <= eps * (Math.Abs(a) + Math.Abs(b) + 1);

    static (Complex, Complex) Pair(double a, double b, double c, double d)
    {
        Complex root = Complex.Sqrt(new Complex((a - d) * (a - d) + 4 * b * c, 0));
        return ((a + d + root) / 2, (a + d - root) / 2);
    }

    public static void Run(double[,] a, double eps, int maxSteps)
    {
        Matrix.Square(a);
        if (eps <= 0 || maxSteps <= 0) throw new ArgumentException("Неверные параметры");
        Matrix.Print("A", a);
        var (q, r) = Decompose(a);
        Matrix.Print("Q", q);
        Matrix.Print("R", r);
        Matrix.Print("Q * R", Matrix.Multiply(q, r));
        Console.WriteLine($"Ошибка QR - A = {Matrix.Difference(Matrix.Multiply(q, r), a):E6}");
        Console.WriteLine($"Ошибка Q^T Q - E = {Matrix.Difference(Matrix.Multiply(Matrix.Transpose(q), q), Matrix.Identity(a.GetLength(0))):E6}");
        var (values, steps) = Eigenvalues(a, eps, maxSteps);
        Console.WriteLine($"epsilon = {eps:E6}; QR-итераций: {steps}");
        for (int i = 0; i < values.Length; i++)
        {
            Complex value = values[i];
            string str = Math.Abs(value.Imaginary) < 1e-9 ? $"{value.Real:F10}" : $"{value.Real:F10} {(value.Imaginary < 0 ? "-" : "+")} {Math.Abs(value.Imaginary):F10}i";
            Console.WriteLine($"lambda{i + 1} = {str}");
        }
    }
}
