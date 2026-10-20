public static class JacobiRotation
{
    public static void Run(double[,] source, double eps, int maxSteps)
    {
        Matrix.Square(source);
        if (eps <= 0 || maxSteps <= 0) throw new ArgumentException("Неверные параметры");
        int n = source.GetLength(0);
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                if (Math.Abs(source[i, j] - source[j, i]) > 1e-10)
                    throw new ArgumentException("Матрица несимметрична");
        double[,] a = (double[,])source.Clone(), v = Matrix.Identity(n);
        int steps = 0;
        double error = OffDiagonal(a);
        Console.WriteLine($"epsilon = {eps:E6}");
        Matrix.Print("A", source);
        while (error > eps && steps < maxSteps)
        {
            int p = 0, q = 1;
            double max = 0;
            for (int i = 0; i < n; i++)
                for (int j = i + 1; j < n; j++)
                    if (Math.Abs(a[i, j]) > max) { max = Math.Abs(a[i, j]); p = i; q = j; }
            if (max == 0) break;
            double app = a[p, p], aqq = a[q, q], apq = a[p, q];
            double tau = (aqq - app) / (2 * apq);
            double t = (tau >= 0 ? 1 : -1) / (Math.Abs(tau) + Math.Sqrt(1 + tau * tau));
            double c = 1 / Math.Sqrt(1 + t * t), s = t * c;
            for (int i = 0; i < n; i++)
            {
                if (i != p && i != q)
                {
                    double aip = a[i, p], aiq = a[i, q];
                    a[i, p] = a[p, i] = c * aip - s * aiq;
                    a[i, q] = a[q, i] = s * aip + c * aiq;
                }
                double vip = v[i, p], viq = v[i, q];
                v[i, p] = c * vip - s * viq;
                v[i, q] = s * vip + c * viq;
            }
            a[p, p] = app - t * apq;
            a[q, q] = aqq + t * apq;
            a[p, q] = a[q, p] = 0;
            steps++;
            error = OffDiagonal(a);
            Console.WriteLine($"Итерация {steps}: внедиагональная норма = {error:E9}");
        }
        if (error > eps) throw new InvalidOperationException("Вращения не сошлись");
        double[] eigenvalues = new double[n];
        double[,] lambda = new double[n, n];
        for (int i = 0; i < n; i++) eigenvalues[i] = lambda[i, i] = a[i, i];
        Console.WriteLine($"Итераций: {steps}; epsilon = {eps:E6}; итоговая норма = {error:E9}");
        Matrix.Print("Собственные значения", eigenvalues);
        for (int j = 0; j < n; j++)
        {
            double[] vector = new double[n];
            for (int i = 0; i < n; i++) vector[i] = v[i, j];
            Matrix.Print($"Собственный вектор {j + 1}", vector);
        }
        Matrix.Print("Матрица V", v);
        Matrix.Print("A * V", Matrix.Multiply(source, v));
        Matrix.Print("V * Lambda", Matrix.Multiply(v, lambda));
        Console.WriteLine($"Ошибка AV - VΛ = {Matrix.Difference(Matrix.Multiply(source, v), Matrix.Multiply(v, lambda)):E6}");
        Console.WriteLine($"Ошибка V^T V - E = {Matrix.Difference(Matrix.Multiply(Matrix.Transpose(v), v), Matrix.Identity(n)):E6}");
    }

    static double OffDiagonal(double[,] a)
    {
        double sum = 0;
        for (int i = 0; i < a.GetLength(0); i++)
            for (int j = i + 1; j < a.GetLength(1); j++) sum += 2 * a[i, j] * a[i, j];
        return Math.Sqrt(sum);
    }
}
