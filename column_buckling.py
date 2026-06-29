def find_critical_load(L, E, A, r, c, e, sigma_allow):
    """
    L: אורך במ"מ
    E: מודול אלסטיות ב-MPa
    A: שטח חתך בממ"ר
    r: רדיוס אינרציה במ"מ
    c: מרחק לסיב קיצוני במ"מ
    e: אקסצנטריות במ"מ
    sigma_allow: מאמץ מותר ב-MPa

    Return: העומס P בניוטון (float)
    """

    def sigma_max(P):
        angle = (L / (2 * r)) * np.sqrt(P / (E * A))
        secant = 1 / np.cos(angle)

        return (P / A) * (1 + (e * c / (r ** 2)) * secant)

    def f(P):
        return sigma_max(P) - sigma_allow

    P_euler = (np.pi ** 2 * E * A * r ** 2) / (L ** 2)

    lower = 0.0
    upper = P_euler * 0.999999

    P_critical = bisect(f, lower, upper, xtol=1e-6, rtol=1e-12)

    return float(P_critical)
