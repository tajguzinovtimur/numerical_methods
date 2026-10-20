public static class Iterations
{
    public static (double[] x, int steps) Jacobi(double[,] a, double[] b, double eps, int maxSteps)
    {
        int n = b.Length;
        double[] x = new double[n], next = new double[n];
        for (int step = 1; step <= maxSteps; step++)
        {
            for (int i = 0; i < n; i++)
            {
                double sum = b[i];
                for (int j = 0; j < n; j++) if (i != j) sum -= a[i, j] * x[j];
                next[i] = sum / a[i, i];
            }
            double delta = Matrix.Difference(x, next);
            Array.Copy(next, x, n);
            if (delta <= eps && Residual(a, x, b) <= eps) return (x, step);
        }
        throw new InvalidOperationException("Простые итерации не сошлись");
    }

    public static (double[] x, int steps) Seidel(double[,] a, double[] b, double eps, int maxSteps)
    {
        int n = b.Length;
        double[] x = new double[n], old = new double[n];
        for (int step = 1; step <= maxSteps; step++)
        {
            Array.Copy(x, old, n);
            for (int i = 0; i < n; i++)
            {
                double sum = b[i];
                for (int j = 0; j < n; j++) if (i != j) sum -= a[i, j] * x[j];
                x[i] = sum / a[i, i];
            }
            if (Matrix.Difference(x, old) <= eps && Residual(a, x, b) <= eps) return (x, step);
        }
        throw new InvalidOperationException("Зейдель не сошёлся");
    }

    public static double Residual(double[,] a, double[] x, double[] b)
    {
        return Matrix.Difference(Matrix.Multiply(a, x), b);
    }

    public static void Run(double[,] a, double[] b, double eps, int maxSteps)
    {
        Matrix.Square(a);
        if (b.Length != a.GetLength(0) || eps <= 0 || maxSteps <= 0)
            throw new ArgumentException("Неверные входные данные");
        int n = b.Length;
        double norm = 0;
        for (int i = 0; i < n; i++)
        {
            if (Math.Abs(a[i, i]) < 1e-13) throw new ArgumentException("Нулевой диагональный элемент");
            double sum = 0;
            for (int j = 0; j < n; j++) if (j != i) sum += Math.Abs(a[i, j] / a[i, i]);
            norm = Math.Max(norm, sum);
        }
        Matrix.Print("A", a);
        Matrix.Print("b", b);
        Console.WriteLine($"epsilon = {eps:E6}; ||B||inf = {norm:F9}");
        if (norm >= 1) Console.WriteLine("Достаточное условие ||B||inf < 1 не выполнено");
        var simple = Jacobi(a, b, eps, maxSteps);
        var seidel = Seidel(a, b, eps, maxSteps);
        var (l, u, p, _) = LU.Decompose(a);
        double[] reference = LU.Solve(l, u, p, b);
        Matrix.Print("Эталонное решение через LU", reference);
        Matrix.Print("Простые итерации", simple.x);
        Console.WriteLine($"Итераций: {simple.steps}; невязка: {Residual(a, simple.x, b):E6}; ошибка: {Matrix.Difference(simple.x, reference):E6}");
        Matrix.Print("Метод Зейделя", seidel.x);
        Console.WriteLine($"Итераций: {seidel.steps}; невязка: {Residual(a, seidel.x, b):E6}; ошибка: {Matrix.Difference(seidel.x, reference):E6}");
        Console.WriteLine(simple.steps == seidel.steps ? "Число итераций одинаковое" : simple.steps < seidel.steps ? "Простые итерации требуют меньше шагов" : "Зейдель требует меньше шагов");
    }
}
